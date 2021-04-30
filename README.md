# SGA_RCT
There are two main modules:
-Design: conducts propensity score estimation and balance check in subgroups; generates the ConnectS plot introduced in Yang et al. (2020)
-Analysis: calculates the weighted subgroup average treatment effect estimates and/or tests of heterogeneous treatment effect across subgroup levels.

Both modules allow users to fit a specific propensity score model, or to provide estimated propensity scores for subsequent analysis. It also offers the option to use flexible machine learning propensity score models (logistic regression, GBM, RFs, restricted cubic spline or pLASSO) and different weighting schemes ( IPW, OW, ATT weights or entropy weights).
Illustrative data examples are provided in both modules. 
