ef train_hmm(train_df):
    # Compute transition counts for states ['c', 'i', 'h', 'M', 'C', 'O', 'n', 'o']
    trans_count = defaultdict(lambda: defaultdict(int))
    obs_count = defaultdict(lambda: defaultdict(int))
    for ann, seq in zip(train_df['ann'], train_df['seq']):
        for i, (pos, next_pos) in enumerate(zip(ann[0:], ann[1:])):
            trans_count[pos][next_pos] += 1
            obs_count[seq[i]][pos] += 1

    # Compute prior counts
    first_count = defaultdict(int)
    for ann in train_df['ann']:
        first_count[ann[0]] += 1

    # Create transition matrix
    size = len(trans_count)
    transition = np.zeros((size, size))
    for i, key in enumerate(trans_count.keys()):
        for j, key2 in enumerate(trans_count.keys()):
            total =  sum(trans_count[key].values())
            transition[i, j] = 1.0*trans_count[key][key2]/total

    # Create prior array
    prior = np.zeros(size)
    for i, key in enumerate(trans_count.keys()):
        prior[i] = first_count[key]
    prior = prior/np.sum(prior)

    # Create observation matrix
    obs_size = len(obs_count)
    observation = np.zeros((size, obs_size))
    for i, key in enumerate(trans_count.keys()):
        for j, key2 in enumerate(obs_count.keys()):
            total = sum(obs_count[key2].values())
            observation[i, j] = 1.0*obs_count[key2][key]/total  
    return trans_count.keys(), transition, prior, obs_count.keys(), observation

