# Signal peptide classifier

# Log

## 30 dec 2017

After further studying previous work in the area, by Nielsen and also the following paper:

[https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2323952/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2323952/)

We decided to also use the annotations in the data during training. We therefore implemented the HMM-package in Biopython. We initializing random probabilities, except for the signal-peptides where probability was set to 1.0 for 'n'. We got the following accuracy using 2-folded cross validation:

```Running k-fold...
Average score  0.912961567445
Deviation +-  0.000376789751319```

Also generated some data files for creating sequence logos.

## 29 dec 2017
* Ran some initial test runs. 

Two parameters were tuned:

1. Number of hidden states
Used 43 hidden states for the positive model and 20 hidden states (chosen quite arbitrarly) for the negative model.
Number of hidden states in the negative model should be investigated further.

Reached an avg. accuracy of about 70% on the validation set using 2-folded cross validation

2. Number of protein sequences 
While Nielsen (1998) truncated all sequences to the first 70 acids, we found the highest accuracy to be yielded using the
first 30 acids. This should be investigated further.

Reached an avg. accuracy of about 85% on the validation set using 2-folded cross validation (30 acids).

## 28 dec 2017

* Wrote code for loading the train/test files.
* Implemented a HMM-model for classifying the signal peptides. 
This was done with inspiration from two articles:

1. [Prediction of signal peptides and signal anchors by a hidden Markov model](https://www.aaai.org/Papers/ISMB/1998/ISMB98-015.pdf)

2. [Prediction of lipoprotein signal peptides in Gram-negative bacteria
](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2323952)

The HMM was implemented using hmmlearn and its [documentation](http://hmmlearn.readthedocs.io/en/latest/tutorial.html).
