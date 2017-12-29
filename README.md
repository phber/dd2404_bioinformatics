# Signal peptide classifier

# Log

## 29 dec 2017
* Ran some initial test runs. 

Two parameters were tuned:

1. Number of hidden states
Used 43 hidden states for the positive model and 20 hidden states (chosen quite arbitrarly) for the negative model.
Reached an accuracy of about 70% on the validation set using 5-folded cross validation. Number of hidden states in the negative model should be investigated further.

2. Number of protein sequences 
While Nielsen (1998) truncated all sequences to the first 70 acids, we found the highest accuracy to be yielded using the
first 30 acids (about 80% accuracy). This should be investigated further.


## 28 dec 2017

1. Wrote code for loading the train/test files.
2. Implemented a HMM-model for classifying the signal peptides. 
This was done with inspiration from two articles:

[Prediction of signal peptides and signal anchors by a hidden Markov model](https://www.aaai.org/Papers/ISMB/1998/ISMB98-015.pdf)

[Prediction of lipoprotein signal peptides in Gram-negative bacteria
](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2323952)

The HMM was implemented using hmmlearn and its [documentation](http://hmmlearn.readthedocs.io/en/latest/tutorial.html).
