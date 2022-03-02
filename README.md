# tsanom

This package contains functions for the detection of anomalies on time series data streams. These methods are appropriate for data with a regular persistant pattern (e.g. daily periodicity or seasonality). THe methods are additionally appropriate for detecting both point (global extreme outliers) and contexual anomalies (local outliers corresponding points that are unusual given their context) and is an online approach.

Two methods are proposed, the first is FDARIMA which combines a short-term model, using an ARIMA model, and long term model, using Functional Data Analysis, to detect both point and contextual anomalies. The second approach, SL Decomposition instead uses a time series decomposition. The long-term structure is first extracted before recovering the short-term structure of the data. Again ARIMA and FDA are used to model the short and long-term components respectively.

*p*-values which are indicators of anomalousness, are used to charactarise whether the newly observed data conforms to the forecasted value. These are calculated using conformal prediction.

This package implements the methods of Elizabeth Riddle-Workman, Marina Evangelou and Niall M. Adams from the work "Combining Forecasts for Improved Anomaly Detection".


## Abstract

Much of the work surrounding anomaly detection in time series data concentrates on finding extreme global outliers known as point anomalies. As time series data often exhibits periodicity or seasonality, it is also important to identify when the series deviates from these behaviours, known as contextual anomalies. This paper presents two anomaly detection procedures, SL Decomposition and FDARIMA, which combine two components; the first capturing the long-term structure of the data such as daily periodic behaviours, the second describing the short-term or local structure of the data. Through modelling these two features separately and using their forecasts, both extreme point anomalies and contextual anomalies can be identified where finding these two types of anomalies is the primary aim of this paper. In addition, the proposed methodology is suitable for data streams where the models continuously adapt to the current behaviour of the series and are computationally fast. The performance of the algorithms is evaluated on both simulated and real data sets and is compared to competitive approaches. The implementation of SL Decomposition and FDARIMA are available in the R package *tsanom*.
