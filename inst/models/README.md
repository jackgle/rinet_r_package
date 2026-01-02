# Place your model and scaler files in this directory

This directory should contain:

1. **rinet_1d.keras** - 1D CNN model for predicting mixture statistics
2. **rinet_2d.keras** - 2D CNN model for predicting mixture statistics  
3. **scaler_1d.pkl** - RobustScaler for 1D model outputs
4. **scaler_2d.pkl** - RobustScaler for 2D model outputs

These files are automatically loaded when you call `predict_rinet_1d()` or 
`predict_rinet_2d()` for the first time.
