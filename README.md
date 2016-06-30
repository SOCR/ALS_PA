# ALS_PA
Amyotrophic Lateral Sclerosis Predictive Analytics (ALS_PA)

In SOCR project uses model-based (e.g., linear models) and model-free (e.g., machine-learning) methods to predict the change of the 
Amyotrophic Lateral Sclerosis Functional Rating Scale (ALSFRS) scores over time. Using independent training data (for model estimation 
and machine learning) and testing data (for assessing the accuracy and reliability of the predictions) we quantify the performance and 
compare the results of different forecasting techniques. The change of the ALSFRS slope between 3 and 12 month of monitoring is 
used as a diagnostic predictive outcome variable. 

Training and testing data of patients and controls from the open Pooled Resource Open-Access ALS Clinical Trials Database (PRO-ACT) archive
are used for demonstration of hte ALS_PA approach.  

Ther R code is released under permissive LGPL license.

Internal validation (using statistical cross-validation), external validation (using the independently provided testing data), 
and separate ALS registry data validation are shown. 
