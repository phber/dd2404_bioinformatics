# Lab notebook

## 2 jan 2018

Final tests were made on the two proteomes. Results seemed to be independent of initialization of the model, so only 1 validation was performed. All squences for the proteom of Bacillus Subtilus was classified by SignalP-TM network only for some uncertain reasons. 

```
Type: PROTEOM of ecoli
Validations : 1
ALL TESTDATA: Average precision, recall, F-score, MCC:
0.654513888889,0.856818181818,0.742125984252,0.716652759202+-0.0,0.0,0.0,0.0
TM TESTDATA: Average precision, recall, F-score, MCC:
0.110344827586,0.761904761905,0.192771084337,0.256856258574+-0.0,0.0,0.0,0.0
NON TM TESTDATA: Average precision, recall, F-score, MCC:
0.837587006961,0.861575178998,0.849411764706,0.827310900299+-0.0,0.0,0.0,0.0


Type: PROTEOM of bacillus
Validations : 1
ALL TESTDATA: Average precision, recall, F-score, MCC:
0.710447761194,0.840989399293,0.770226537217,0.755222944136+-0.0,0.0,0.0,0.0
TM TESTDATA: Average precision, recall, F-score, MCC:
+-
NON TM TESTDATA: Average precision, recall, F-score, MCC:
+-
```

## 1 jan 2018

1. To balance the precision and recall for TM-testdata, we decided to add a treshhold for scores between the models. If the score diferrence between the model is less than 2, the prediction is 0. This made the precision and recall for TM-data more even (see res file).

2. We disallow some transitions in the model that are impossible (e.g, n -> c).

3. Add functions for reading proteomes. Also added functions for reading SignalP 4.1 output files, to use as a baseline error.

Fasta sequences were submitted to SignalP 4.1 and outputs saved in res/proteomes. Some test runs on Ecoli were made to measure the F-score of the classifier (see res file).

4. Added metrics for MCC.

## 31 dec 2017

Changed to performance measures to precision, recall and F-score. Also changed testing to test only TM or non-TM proteins, or both. The results were (using 5-folded validation):
```
Type: CROSS VALIDATION
Validations: 5
ALL TESTDATA: Average precision, recall, F-score:
[0.89952194928962004, 0.92475784707141406, 0.91128394619236386]
TM TESTDATA: Average precision, recall, F-score:
[0.5696525105388679, 0.89361184426298856, 0.69494933638859258]
NON TM TESTDATA: Average precision, recall, F-score:
[0.91996094358046965, 0.93647088476833673, 0.92801386388
```
These can also be found in the "res" folder. 

## 30 dec 2017

After further studying previous work in the area, by Nielsen and also the following paper:

[https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2323952/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2323952/)

We decided to also use the annotations in the data during training. We therefore implemented the HMM-package in Biopython. We initializing random probabilities, except for the signal-peptides where probability was set to 1.0 for 'n'. We got the following accuracy using 2-folded cross validation:

```
Running k-fold...
Average score  0.912961567445
Deviation +-  0.000376789751319
```

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
