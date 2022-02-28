# tsanom

This package contains functions for the detection of anomalies on time series data streams. These methods are appropriate for data with a regular persistant pattern (e.g. daily periodicity or seasonality). THe methods are additionally appropriate for detecting both point (global extreme outliers) and contexual anomalies (local outliers corresponding points that are unusual given their context) and is an online approach.

Two methods are proposed, the first is FDARIMA which combines a short-term model, using an ARIMA model, and long term model, using Functional Data Analysis, to detect both point and contextual anomalies. The second approach, SL Decomposition instead uses a time series decomposition. The long-term structure is first extracted before recovering the short-term structure of the data. Again ARIMA and FDA are used to model the short and long-term components respectively.

*p*-values which are indicators of anomalousness, are used to charactarise whether the newly observed data conforms to the forecasted value. These are calculated using conformal prediction.

This package implements the methods of Elizabeth Riddle-Workman, Marina Evangelou and Niall M. Adams from the work "Combining Forecasts for Improved Anomaly Detection".
