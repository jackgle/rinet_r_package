"""Convert a trained RINet Keras model to TorchScript for the R package.

Usage:
    python data-raw/convert_keras_to_torchscript.py MODEL.keras OUT.pt [--scaler SCALER.pkl]

Builds the equivalent torch network from the Keras layer configs (Conv1D/2D and
MaxPooling1D/2D with 'same' padding, Flatten, Dense, Dropout), copies the
weights, checks the outputs against Keras on random inputs, and saves a
TorchScript trace. With --scaler, also prints the RobustScaler constants in the
form used by .scaler_params in R/model_functions.R.

Requires keras (TensorFlow backend), torch, numpy, and scikit-learn (--scaler).
"""
import argparse
import os
import pickle

os.environ.setdefault('TF_CPP_MIN_LOG_LEVEL', '3')

import keras
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F


class SamePadMaxPool(nn.Module):
    """Keras MaxPooling with padding='same': pad with -inf, then pool."""

    def __init__(self, ndim, pool_size, pads):
        super().__init__()
        self.pool = (nn.MaxPool1d if ndim == 1 else nn.MaxPool2d)(pool_size)
        self.pads = pads  # F.pad order: last spatial dim first

    def forward(self, x):
        return self.pool(F.pad(x, self.pads, value=float('-inf')))


class ChannelsLastFlatten(nn.Module):
    """Flatten in Keras's channels-last element order."""

    def __init__(self, ndim):
        super().__init__()
        self.order = [0, 2, 1] if ndim == 1 else [0, 2, 3, 1]

    def forward(self, x):
        return torch.flatten(x.permute(self.order), start_dim=1)


def activation(name):
    if name == 'relu':
        return nn.ReLU()
    if name == 'linear':
        return nn.Identity()
    raise ValueError(f"Unsupported activation: {name}")


def build_torch_model(keras_model):
    in_shape = keras_model.input_shape[1:]  # channels-last, e.g. (100, 1)
    ndim = len(in_shape) - 1
    spatial, channels = list(in_shape[:-1]), in_shape[-1]
    layers = []

    for layer in keras_model.layers:
        kind, cfg = type(layer).__name__, layer.get_config()

        if kind in ('InputLayer', 'Dropout'):  # dropout is inactive at inference
            continue

        if kind in ('Conv1D', 'Conv2D'):
            k = cfg['kernel_size']
            assert cfg['padding'] == 'same' and all(s == 1 for s in cfg['strides'])
            assert all(i % 2 == 1 for i in k), "even kernels need asymmetric padding"
            conv = (nn.Conv1d if ndim == 1 else nn.Conv2d)(
                channels, cfg['filters'], k, padding=[i // 2 for i in k])
            w, b = layer.get_weights()
            # Keras (k..., in, out) -> torch (out, in, k...)
            conv.weight.data = torch.from_numpy(np.moveaxis(w, [-1, -2], [0, 1]).copy())
            conv.bias.data = torch.from_numpy(b.copy())
            layers += [conv, activation(cfg['activation'])]
            channels = cfg['filters']

        elif kind in ('MaxPooling1D', 'MaxPooling2D'):
            p = cfg['pool_size']
            assert cfg['padding'] == 'same' and tuple(cfg['strides']) == tuple(p)
            out = [-(-s // q) for s, q in zip(spatial, p)]
            pads = []
            for o, q, s in reversed(list(zip(out, p, spatial))):
                need = (o - 1) * q + q - s
                pads += [need // 2, need - need // 2]
            layers.append(SamePadMaxPool(ndim, p, pads))
            spatial = out

        elif kind == 'Flatten':
            layers.append(ChannelsLastFlatten(ndim))
            channels, spatial = int(np.prod(spatial)) * channels, []

        elif kind == 'Dense':
            w, b = layer.get_weights()
            lin = nn.Linear(w.shape[0], w.shape[1])
            lin.weight.data = torch.from_numpy(w.T.copy())
            lin.bias.data = torch.from_numpy(b.copy())
            layers += [lin, activation(cfg['activation'])]

        else:
            raise ValueError(f"Unsupported layer type: {kind}")

    return nn.Sequential(*layers).eval(), in_shape


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('model')
    parser.add_argument('out')
    parser.add_argument('--scaler')
    args = parser.parse_args()

    keras_model = keras.models.load_model(args.model, compile=False)
    torch_model, in_shape = build_torch_model(keras_model)

    x = np.random.default_rng(0).uniform(size=(16, *in_shape)).astype('float32')
    x_torch = torch.from_numpy(np.moveaxis(x, -1, 1).copy())  # channels-first
    with torch.no_grad():
        traced = torch.jit.trace(torch_model, x_torch[:1])
        max_diff = np.abs(keras_model.predict(x, verbose=0) - traced(x_torch).numpy()).max()
    print(f"max |keras - torch| on random inputs: {max_diff:.2e}")
    assert max_diff < 1e-4, "TorchScript output does not match Keras"
    traced.save(args.out)
    print(f"saved {args.out}")

    if args.scaler:
        with open(args.scaler, 'rb') as f:
            scaler = pickle.load(f)
        fmt = lambda v: ', '.join(repr(float(i)) for i in v)
        print(f"center = c({fmt(scaler.center_)}),")
        print(f"scale  = c({fmt(scaler.scale_)})")


if __name__ == '__main__':
    main()
