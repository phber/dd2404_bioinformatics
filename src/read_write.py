import os
from Bio import SeqIO
import pandas as pd
import datetime
from Bio import Alphabet
import numpy as np

""" Paths for training, test data, proteomes """
BASE_DIR = os.path.join(os.path.dirname(__file__ ), '..')
DATA_DIR = os.path.join(BASE_DIR, 'data')
RES_DIR = os.path.join(BASE_DIR, 'res')
PROTEOM_DATA_DIR = os.path.join(DATA_DIR, 'proteomes')
PROTEOM_RES_DIR = os.path.join(RES_DIR, 'proteomes')
NEGATIVE_TM_DIR = os.path.join(DATA_DIR, 'negative_examples', 'tm')
POSITIVE_TM_DIR = os.path.join(DATA_DIR, 'positive_examples', 'tm')
NEGATIVE_NON_TM_DIR = os.path.join(DATA_DIR, 'negative_examples', 'non_tm')
POSITIVE_NON_TM_DIR = os.path.join(DATA_DIR, 'positive_examples', 'non_tm')

#Index for truncation of sequences
TRUNC_IDX = 70

"""Alphabet for labels and sequences"""
seq_alphabet = Alphabet.SingleLetterAlphabet()
seq_alphabet.letters = ['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y', 'X']
state_alphabet_pos = Alphabet.SingleLetterAlphabet()
state_alphabet_pos.letters = ['i', 'c', 'h', 'C', 'O', 'n', 'M', 'O', 'o']
state_alphabet_neg = Alphabet.SingleLetterAlphabet()
state_alphabet_neg.letters = ['i', 'M', 'O', 'o']

""" Loads all fasta sequences from a folder into a DataFrame"""
def load_dir(folder):
    data = {'descr': [], 'id': [], 'seq': [], 'ann' : []}
    for _f in os.listdir(folder):
        f = os.path.join(folder, _f)
        for record in SeqIO.parse(f, 'fasta'):
            seq, ann = record.seq.split('#')
            seq.alphabet = seq_alphabet
            data['seq'].append(seq[:TRUNC_IDX])
            data['ann'].append(ann[:TRUNC_IDX])
            data['descr'].append(record.description)
            data['id'].append(record.id)
    return pd.DataFrame(data)

"""Returns sequences from a df"""
def get_secrecords(df, lenlimit = 1000):
    seqs = []
    from collections import defaultdict
    inits = defaultdict(int)
    for i, row in df.iterrows():
        inits[str(row['ann'])[0]] += 1
        if len(row['seq']) < lenlimit:
            seqs.append(SeqIO.SeqRecord(row['seq'], id=row['id']))
    print inits
    return seqs

""" Fills a dataframe with data using FOLDER CONSTANTS"""
def fill_train_df(write = False):
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
    res = pd.concat([neg_tm, pos_tm, neg_non_tm, pos_non_tm]).reset_index(drop=True)
    if write:
        seqs_pos = get_secrecords(pos_tm) + get_secrecords(pos_non_tm)
        with open("res/logo_positive.fasta", "w") as output_handle:
            SeqIO.write(seqs_pos, output_handle, "fasta")
    return res

"""Writes scores to res folder"""
def store_result(scores_all, scores_tm, scores_non_tm, name, validations):
    with open(os.path.join(RES_DIR, str(datetime.date.today()) + '.txt'), 'a') as f:
        f.write('\n')
        f.write('Type: ' + name)
        f.write('\n')
        f.write('Validations : ' + str(validations))
        f.write('\n')
        f.write('ALL TESTDATA: Average precision, recall, F-score, MCC, (Cleave Acc.):\n')
        f.write(','.join([str(np.mean(y)) for y in zip(*scores_all)]))
        f.write('+-')
        f.write(','.join([str(np.std(y)) for y in zip(*scores_all)]))
        f.write('\n')
        f.write( 'TM TESTDATA: Average precision, recall, F-score, MCC, (Cleave Acc.):\n')
        f.write(','.join([str(np.mean(y)) for y in zip(*scores_tm)]))
        f.write('+-')
        f.write(','.join([str(np.std(y)) for y in zip(*scores_tm)]))
        f.write('\n')
        f.write( 'NON TM TESTDATA: Average precision, recall, F-score, MCC (Cleave Acc.):\n')
        f.write(','.join([str(np.mean(y)) for y in zip(*scores_non_tm)]))
        f.write('+-')
        f.write(','.join([str(np.std(y)) for y in zip(*scores_non_tm)]))
        f.write('\n')


"""Splits a proteom into fasta files of max length 2000 for SignalP compatability"""
def split_file(filename):
    f = os.path.join(PROTEOM_DATA_DIR, filename)
    records = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))
    name_index = 1
    ix = range(0, len(records), 2000)
    ix.append(len(records))
    for i in range(1, len(ix)):
        new_f = os.path.join(PROTEOM_RES_DIR, str(name_index) + '_' + filename)
        with open(new_f, 'w') as output_handle:
           SeqIO.write(records.values()[ix[i-1]:ix[i]], output_handle, 'fasta')
        name_index += 1

"""Loads a Signal P result file, named IDX_NAME_signalp.fasta"""
def load_signalp(proteom, stopidx, startidx = 1):
    data = {'descr': [], 'id': [], 'label' : [], 'seq' : [], 'tm' : []}
    for i in range(startidx, stopidx + 1):
        #Load signal P labels
        f = os.path.join(PROTEOM_RES_DIR, str(i) + '_' + proteom + '_signalp.fasta')
        for r in SeqIO.parse(f, 'fasta'):
            data['descr'].append(r.description)
            data['id'].append(r.id)
            if "SP='YES'" in str(r.seq):
                data['label'].append(1)
            elif "SP='NO'" in str(r.seq):
                data['label'].append(0)
            else:
                raise ValueError("Sequence missing SP-label, e.g SP='YES'")
            if 'Networks=SignalP-TM' in str(r.seq):
                data['tm'].append(True)
            elif 'Networks=SignalP-noTM' in str(r.seq):
                data['tm'].append(False)
            else: 
                raise ValueError("Sequence missing TM-label, e.g Networks=SignalP-noTM")
        # Load data sequences 
        seq_f = os.path.join(PROTEOM_DATA_DIR, str(i) + '_' + proteom + '.fasta')
        for ix, r in enumerate(SeqIO.parse(seq_f, 'fasta')):
            signalp_id = r.id.replace('|', '_')
            assert data['id'][2000*(i-1) + ix] == signalp_id, 'ID mismatch between SignalP and sequence files.'
            data['seq'].append(r.seq[:TRUNC_IDX])
    return pd.DataFrame(data)

split_file('bacillus.fasta')
