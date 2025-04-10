# parTimeROC 0.1.1 (2025-04-10) 

Small improvements and fixes in `timeroc_gof` and `timeroc_predict` method. Added acknowledgement section in the DESCRIPTION file.

# parTimeROC 0.1.0 (2024-06-13) 

This is a major release of the `parTimeROC` package which enables running the Time-Dependent Receiver Operating Characteristics (ROC) using parametric approaches. Two models were used which are based on the Proportional Hazard and copula functions.

## Features

Several methods prepared in this package are as follows:
1. `timeroc_obj()` - To create a `TimeROC` object.
2. `rtimeroc()` - To simulate random data based on the chosen model.
3. `timeroc_fit()` - To estimate model's parameter using either the frequentist or Bayesian.
4. `timeroc_gof()` - To check the model's goodness-of-fit.
5. `timeroc_predict()` - To calculate time-dependent ROC curve at selected time point.
6. `timeroc_auc()` - To calculate the area under the time-dependent ROC curve.
7. `rate_change()` - To calculate the rate change of the time-dependent ROC curve.
