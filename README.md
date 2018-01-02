# Signal peptide discrimination using a HMM with known states

A classifier that is able to detect signal peptide sequences in a proteom, using two HMMs that are trained with annotated sequences (known states).

## Instructions

1. Load training data:

```df = fill_train_df()```

2. To run cross validation on training data use:

```cross_validate(df, validations = 5, write = True)```

3. To run a performance test on a proteome (labelled by SignalP 4.1) use:

```test_proteom(df, 'ecoli', stopidx = 3, write = True)```

All results are stored in the 'res' folder.
