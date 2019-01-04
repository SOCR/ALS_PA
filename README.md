# ALS_PA
Amyotrophic Lateral Sclerosis Predictive Analytics (**ALS_PA**)

This [SOCR](http://socr.umich.edu/) project uses model-based (e.g., linear models) and model-free (e.g., machine-learning) methods to predict the change of the Amyotrophic Lateral Sclerosis Functional Rating Scale (ALSFRS) scores over time. Using independent training data (for model estimation and machine learning) and independent testing data (for assessing the accuracy and reliability of the predictions) we quantify the algorithmic performance and compare the results of different forecasting techniques. By training models of the change of the ALSFRS slope from baseline to 3 (learning/training) and subsequently testing the forecasts on the ALSFRS change between 3 and 12 months allows us to quantify the performance of the predictors to accurately diagnose the disease severity and progression over time. 

Training and testing data of patients and controls from the open Pooled Resource Open-Access ALS Clinical Trials Database (PRO-ACT) archive are used for demonstration of the ALS_PA approach.  

In support of [open science](https://en.wikipedia.org/wiki/Open_science), the complete **R code** used in this study is shared here, under [permissive LGPL license](http://socr.umich.edu/html/SOCR_CitingLicense.html), to enable independent result reproducibility, community validation, and transdisciplinary team science.

Internal validation (using statistical cross-validation), external validation (using the independently provided testing data), 
and separate ALS registry data validation are also demonstrated. 

## References

* [SOCR Website](http://www.socr.umich.edu).
* Tang, M., Gao, C, Goutman, SA, Kalinin, A, Mukherjee, B, Guan, Y, and Dinov, ID. (2018) [Model-Based and Model-Free Techniques for Amyotrophic Lateral Sclerosis Diagnostic Prediction and Patient Clustering](https://doi.org/10.1007/s12021-018-9406-9), Neuroinformatics, 1-15, [DOI: 10.1007/s12021-018-9406-9](https://doi.org/10.1007/s12021-018-9406-9).
* Huang Z, Zhang H, Boss J, Goutman SA, Mukherjee B, Dinov ID, Guan, Y. (2017) [Complete hazard ranking to analyze right-censored data: An ALS survival study](https://doi.org/10.1371/journal.pcbi.1005887). PLoS Comput Biol 13(12): e1005887, [DOI: 10.1371/journal.pcbi.1005887](https://doi.org/10.1371/journal.pcbi.1005887).
