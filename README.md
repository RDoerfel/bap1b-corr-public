# Correlation of Rates of Change in Longitudinal Studies with Dropout

Tutorial accompanying the paper *Estimating the correlation between rates of change in longitudinal studies in presence of dropout*.

The paper shows that the common slopes-first estimator (fitting per-subject OLS regressions and correlating the resulting slopes) is biased under independent dropout whenever the cross-covariance between two longitudinal signals is not constant over time. We propose an unbiased alternative that estimates the joint covariance matrix before applying the contrast.

We provide a tutorial for comparing these estimators in the form of an R script (`tutorial_script.R`) and a Jupyter notebook for Google Colab (`tutorial_colab.ipynb`).
