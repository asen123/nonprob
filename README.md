# nonprob

This repository contains codes for the case study in the paper "Improving measurement error and representativeness in nonprobability surveys". 
The paper is available on arXiv at the following link: 
https://arxiv.org/abs/2410.18282.

The original PEW survey data and detailed report are available at the following link: 
https://www.pewresearch.org/methods/2023/09/07/comparing-two-types-of-online-survey-samples/

The files "code step 1.R", "code step 2.R" and "code step 3.R" are the R codes, which need to be run sequentially. Following is a summary of the codes:

code step 1 : This code analyzes the 6 Pew surveys to create IPW estimates, bias-corrected estimates and composite estimates for the case when bias is considered known.
code step 2 : This code creates the aforementioned estimates for the case when bias is considered unknonwn, using modeling on probability sample.
code step 3 : This code performs sampling from probability survey for the purpose of comparison of composite estimators to those only from probability surveys.
