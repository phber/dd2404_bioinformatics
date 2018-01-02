import pandas as pd
from sklearn.model_selection import KFold
import numpy as np
from Bio.HMM import MarkovModel, Trainer
from sklearn.metrics import confusion_matrix
from sklearn.utils import shuffle
from sklearn.metrics import matthews_corrcoef
from read_write import seq_alphabet, state_alphabet_pos, state_alphabet_neg
from read_write import load_dir, load_signalp, fill_train_df, store_result
from Bio.Seq import Seq

"""Train a HMM using ML-estimation"""
def fit_model(train_df, positive = True):
    training_seqs = []
    if positive:
        builder = MarkovModel.MarkovModelBuilder(state_alphabet_pos, seq_alphabet)
        builder.allow_transition('n', 'n')
        builder.allow_transition('n', 'h')
        builder.allow_transition('h', 'h')
        builder.allow_transition('h', 'c')
        builder.allow_transition('c', 'c')
        builder.allow_transition('c', 'C')
        builder.allow_transition('C', 'O')
        builder.allow_transition('O', 'O')
        builder.allow_transition('C', 'o')
        builder.allow_transition('o', 'o')
        builder.allow_transition('o', 'M')
        builder.allow_transition('M', 'M')
        builder.allow_transition('M', 'i')
        builder.allow_transition('i', 'i')
        builder.set_random_probabilities()
        builder.set_initial_probabilities({'n' : 1})
    else: 
        builder = MarkovModel.MarkovModelBuilder(state_alphabet_neg, seq_alphabet)
        builder.allow_all_transitions()
        builder.destroy_transition('O', 'i')
        builder.destroy_transition('O', 'M')
        builder.destroy_transition('O', 'o')
        builder.set_initial_probabilities({'o' : 0.1, 'i' : 0.4, 'O': 0.5})
    training_seqs = []
    for i, row in train_df.iterrows():
        ann = row['ann']            
        ann.alphabet = state_alphabet_pos if positive else state_alphabet_neg
        train_seq = Trainer.TrainingSequence(row['seq'], ann)
        training_seqs.append(train_seq)
    model = builder.get_markov_model()
    trainer = Trainer.KnownStateTrainer(model)
    return trainer.train(training_seqs)

"""Prediction scores"""
def perf_measure(y_true, y_pred):
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    rec = 1.0*tp/(tp+fn)
    prec = 1.0*tp/(tp+fp)
    f = 2*prec*rec/(prec+rec)
    mcc = matthews_corrcoef(y_true, y_pred)
    return (prec, rec, f, mcc)

"""Computes test scores for a df with known labels"""
def hmm_test(pos_model, neg_model, test_df):
    predictions = []

    for i, row in test_df.iterrows():
        seq = str(row['seq'])
        seq = seq.replace('U', '')
        pos_seq, pos_score = pos_model.viterbi(seq, state_alphabet_pos)
        neg_sec, neg_score = neg_model.viterbi(seq, state_alphabet_neg)
        if abs(pos_score - neg_score) < 2:
            predictions.append(0)
        elif pos_score > neg_score:
            predictions.append(1)
        else:
            predictions.append(0)
    print np.sum(predictions)
    return perf_measure(list(test_df['label']), predictions)

"""Helper function for splitting dataframe by binary group"""
def split_df(df, group='label'):
     #Split training into negative and positive
    gb = df.groupby(group)
    groups = [gb.get_group(x) for x in gb.groups]
    train_df_neg = groups[0]
    train_df_pos = groups[1]
    return train_df_neg, train_df_pos

"""Get cross validation scores for training data"""
def cross_validate(df, write = True, validations = 5):
    if len(df) == 0:
        raise ValueError('No training data was found.')
    kf = KFold(n_splits = validations, shuffle = True, random_state = 13)
    scores_all = []
    scores_non_tm = []
    scores_tm = []
    for train_index, test_index in kf.split(df):
        print 'Running k-fold...'
        # Split into train and test
        train_df = df.iloc[train_index]
        test_df =  df.iloc[test_index]

        #Split training into negative and positive
        train_df_neg, train_df_pos = split_df(train_df, 'label')

        #Train models
        pos_model = fit_model(train_df_pos, True)
        neg_model = fit_model(train_df_neg, False)

        #Run tests
        test_df_non_tm, test_df_tm = split_df(test_df, 'tm')
        scores_all.append(hmm_test(pos_model, neg_model, test_df))
        scores_tm.append(hmm_test(pos_model, neg_model, test_df_tm))
        scores_non_tm.append(hmm_test(pos_model, neg_model, test_df_non_tm))
    if write:
        store_result(scores_all, scores_tm, scores_non_tm, 'CROSS VALIDATION', validations)
    else: 
        print 'ALL TESTDATA', scores_all


"""Stats for proteom df"""
def stats(df):
    test_df_non_tm, test_df_tm = split_df(df, 'tm')
    print 'Non-TM count, TM count:'
    print len(test_df_non_tm), len(test_df_tm)

    print 'Negative TM count, Positive TM count:'
    t1, t2 = split_df(test_df_tm)
    print len(t1), len(t2)

    print 'Negative Non-TM count, Positive Non-TM count:'
    t1, t2 = split_df(test_df_non_tm)
    print len(t1), len(t2)

"""Test a model HMM on a proteom with signalP resulsts"""
def test_proteom(train_df, proteom, stopidx, startidx = 1, write = True):
    #Load proteom
    prot_df = load_signalp(proteom, stopidx, startidx)
    if proteom != 'bacillus':
        test_df_non_tm, test_df_tm = split_df(prot_df, 'tm')
        stats(prot_df)
    else:
        print len(prot_df), np.sum(prot_df['label'])
    scores_all = []
    scores_non_tm = []
    scores_tm = []

    train_df_neg, train_df_pos = split_df(train_df)

    pos_model = fit_model(train_df_pos, True)
    neg_model = fit_model(train_df_neg, False)
    scores_all.append(hmm_test(pos_model, neg_model, prot_df))
    if proteom != 'bacillus':
        scores_tm.append(hmm_test(pos_model, neg_model, test_df_tm))
        scores_non_tm.append(hmm_test(pos_model, neg_model, test_df_non_tm))

    if write:
        print scores_all
        store_result(scores_all, scores_tm, scores_non_tm, 'PROTEOM of ' + proteom, 1)
    else: 
        print 'ALL TESTDATA', np.mean(scores_all)
        print 'TM TESTDATA', np.mean(scores_tm)
        print 'NON TM TESTDATA', np.mean(scores_non_tm)

df =  fill_train_df()
#test_proteom(df, 'bacillus', stopidx = 3, write = True)
#cross_validate(df, validations = 2, write = True)