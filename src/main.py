from Bio import SeqIO
import os
import pandas as pd
from sklearn.model_selection import KFold
import numpy as np
from collections import defaultdict
from hmmlearn import hmm

BASE_DIR = os.path.join(os.path.dirname(__file__ ), '..')
DATA_DIR = os.path.join(BASE_DIR, 'data')
NEGATIVE_TM_DIR = os.path.join(DATA_DIR, 'negative_examples', 'tm')
POSITIVE_TM_DIR = os.path.join(DATA_DIR, 'positive_examples', 'tm')
NEGATIVE_NON_TM_DIR = os.path.join(DATA_DIR, 'negative_examples', 'non_tm')
POSITIVE_NON_TM_DIR = os.path.join(DATA_DIR, 'positive_examples', 'non_tm')

""" Loads all fasta sequences from a folder into a DataFrame"""
def load_dir(folder):
    data = {'descr': [], 'id': [], 'seq': [], 'ann' : []}
    for _f in os.listdir(folder):
        f = os.path.join(folder, _f)
        for record in SeqIO.parse(f, 'fasta'):
            data['descr'].append(record.description)
            data['id'].append(record.id)
            seq, ann = record.seq.split('#')
            data['seq'].append(str(seq)[:30])
            data['ann'].append(str(ann))
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

"""Converts a Sequence to numeric format"""
def get_binary(df, obs_keys):
    if type(df) == str:
        res = []
        for s in df:
            if s in obs_keys:
                res.append(obs_keys.index(s)) 
        return res

    binseq = []
    for seq in df['seq']:
        res = []
        for s in seq:
            if s in obs_keys:
                res.append(obs_keys.index(s))
        binseq.append(res)
    return binseq

"""Format df for hmm"""
def getX(df):
    binseqs = np.array(df['binseq'])
    X = np.array([binseqs[0]]).T
    lens = [len(X)]
    for binseq in binseqs[1:]:
        lens.append(len(binseq))
        X = np.concatenate([X, np.array([binseq]).T])
    return X, lens

"""Remove sequence after cleave site"""
def remove_00s(train_df):
    for i, row in train_df.iterrows():
        if 'C' in row['ann']:
            cleave_index = row['ann'].index('C')
            train_df.loc[i, 'seq'] = row['seq'][:cleave_index]
            train_df.loc[i, 'ann'] = row['ann'][:cleave_index]

"""Train a HMM using baum-welch"""
def fit_model(train_df, states):
    model = hmm.MultinomialHMM(n_components=states)
    X, lengths = getX(train_df)
    model.fit(X, lengths)
    return model

"""Get cross validation scores"""
def run_hmm(df):
    kf = KFold(n_splits = 2, shuffle = True, random_state = 2)
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

        #Compute number of states from training data
        trans_keys_pos, obs_keys_pos = calc_states(train_df_pos)
        trans_keys_neg, obs_keys_neg = calc_states(train_df_neg)

        #Convert sequences to integers
        train_df_pos['binseq'] = get_binary(train_df_pos, obs_keys_pos)
        train_df_neg['binseq'] = get_binary(train_df_neg, obs_keys_neg)

        pos_model = fit_model(train_df_pos, 43)
        neg_model = fit_model(train_df_neg, 20) 

        predictions = []
        #Compute test scores
        for seq in test_df['seq']:
            seq_pos = get_binary(seq, obs_keys_pos)
            seq_neg = get_binary(seq, obs_keys_neg)
            pos_score = pos_model.score([seq_pos])
            neg_score = neg_model.score([seq_neg])
            if pos_score > neg_score:
                predictions.append(1)
            else:
                predictions.append(0)
        score = 1.0*np.sum(np.array(test_df['label']) == predictions)/len(predictions)
        #Predict peptide or non peptide
        scores.append(score)
    print 'Average score ', np.mean(scores)
    print 'Deviation +- ', np.std(scores)

df =  fill_df()
run_hmm(df)