from Bio import SeqIO
from Bio import Alphabet
import os
import pandas as pd
from sklearn.model_selection import KFold
import numpy as np
from collections import defaultdict
from Bio.HMM import MarkovModel, Trainer
from sklearn.metrics import confusion_matrix

""" Paths for training and test data """
BASE_DIR = os.path.join(os.path.dirname(__file__ ), '..')
DATA_DIR = os.path.join(BASE_DIR, 'data')
NEGATIVE_TM_DIR = os.path.join(DATA_DIR, 'negative_examples', 'tm')
POSITIVE_TM_DIR = os.path.join(DATA_DIR, 'positive_examples', 'tm')
NEGATIVE_NON_TM_DIR = os.path.join(DATA_DIR, 'negative_examples', 'non_tm')
POSITIVE_NON_TM_DIR = os.path.join(DATA_DIR, 'positive_examples', 'non_tm')

"""Alphabet for labels and sequences"""
seq_alphabet = Alphabet.SingleLetterAlphabet()
seq_alphabet.letters = ['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y', 'X']
state_alphabet_pos = Alphabet.SingleLetterAlphabet()
state_alphabet_pos.letters = ['c', 'h', 'C', 'O', 'n','i', 'M', 'O', 'o']
state_alphabet_neg = Alphabet.SingleLetterAlphabet()
state_alphabet_neg.letters = ['i', 'M', 'O', 'o']

""" Loads all fasta sequences from a folder into a DataFrame"""
def load_dir(folder):
    data = {'descr': [], 'id': [], 'seq': [], 'ann' : []}
    for _f in os.listdir(folder):
        f = os.path.join(folder, _f)
        for record in SeqIO.parse(f, 'fasta'):
            data['descr'].append(record.description)
            data['id'].append(record.id)
            seq, ann = record.seq.split('#')
            seq.alphabet = seq_alphabet
            data['seq'].append(seq[:70])
            data['ann'].append(ann[:70])
    return pd.DataFrame(data)

""" Fills a dataframe with data using FOLDER CONSTANTS"""
def fill_df():
    neg_tm = load_dir(NEGATIVE_TM_DIR)
    neg_tm['label'] = 0
    neg_tm['tm'] = True

    pos_tm = load_dir(POSITIVE_TM_DIR)
    pos_tm['label'] = 1
    pos_tm['tm'] = True

    neg_non_tm = load_dir(NEGATIVE_NON_TM_DIR)
    neg_non_tm['label'] = 0
    neg_non_tm['tm'] = False

    pos_non_tm = load_dir(POSITIVE_NON_TM_DIR)
    pos_non_tm['label'] = 1
    pos_non_tm['tm'] = False
    return pd.concat([neg_tm, pos_tm, neg_non_tm, pos_non_tm]).reset_index(drop=True)

"""Get possible transition states and emissions"""
def calc_states(train_df):
    trans_count = defaultdict(int)
    obs_count = defaultdict(int)
    for ann, seq in zip(train_df['ann'], train_df['seq']):
        for a, s in zip(ann, seq):
            trans_count[a] += 1
            if s != 'X':
                obs_count[s] += 1
    return trans_count.keys(), obs_count.keys()

"""Train a HMM using baum-welch"""
def fit_model(train_df, positive = True):
    training_seqs = []
    if positive:
        builder = MarkovModel.MarkovModelBuilder(state_alphabet_pos, seq_alphabet)
        builder.allow_all_transitions()
        builder.set_random_probabilities()
        builder.set_initial_probabilities({'n' : 1})
    else: 
        builder = MarkovModel.MarkovModelBuilder(state_alphabet_neg, seq_alphabet)
        builder.allow_all_transitions()
        builder.set_random_probabilities()
    training_seqs = []
    for i, row in train_df.iterrows():
        ann = row['ann']
        ann.alphabet = state_alphabet_pos if positive else state_alphabet_neg
        train_seq = Trainer.TrainingSequence(row['seq'], ann)
        training_seqs.append(train_seq)
    model = builder.get_markov_model()
    trainer = Trainer.KnownStateTrainer(model)
    trainer.train(training_seqs)
    return trainer.train(training_seqs)

"""Prediction scores"""
def perf_measure(y_true, y_pred):
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    print tn, fp, fn, tp

"""Get cross validation scores"""
def run_hmm(df):
    kf = KFold(n_splits = 2, shuffle = True, random_state = 4)
    scores = []
    for train_index, test_index in kf.split(df):
        print 'Running k-fold...'
        # Split into train and test
        train_df = df.iloc[train_index]
        test_df =  df.iloc[test_index]

        #Split training into negative and positive
        gb = train_df.groupby('label')
        groups = [gb.get_group(x) for x in gb.groups]
        train_df_neg = groups[0]
        train_df_pos = groups[1]
        pos_model = fit_model(train_df_pos, True)
        neg_model = fit_model(train_df_neg, False) 
        predictions = []
        #Compute test scores
        for i, row in test_df.iterrows():
            seq = row['seq']
            pos_seq, pos_score = pos_model.viterbi(seq, state_alphabet_pos)
            neg_sec, neg_score = neg_model.viterbi(seq, state_alphabet_neg)
            if pos_score > neg_score:
                predictions.append(1)
            else:
                predictions.append(0)
        score = 1.0*np.sum(np.array(test_df['label']) == predictions)/len(predictions)
        print score
        scores.append(score)
        perf_measure(list(test_df['label']), predictions)
    print 'Average score ', np.mean(scores)
    print 'Deviation +- ', np.std(scores)

df =  fill_df()
run_hmm(df)