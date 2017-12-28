# Signal peptide classifier

# Log

## 28 dec 2017

1. Wrote code for loading the train/test files.
2. Implemented a HMM-model for classifying the signal peptides. 
This was done with inspiration from two articles:

[Prediction of signal peptides and signal anchors by a hidden Markov model][https://www.aaai.org/Papers/ISMB/1998/ISMB98-015.pdf]
[Prediction of lipoprotein signal peptides in Gram-negative bacteria
][https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2323952]

The HMM was implemented using hmmlearn and its [documentation](http://hmmlearn.readthedocs.io/en/latest/tutorial.html).
