import os
import numpy   as np
import pandas  as pd
import sklearn as sk
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA, SparsePCA



def import_data(dirstr, round):
    # Read subject data
    sids, data = read_response_data(dirstr)
    
    # Get subject attention check performance
    sids, data = quality_check(sids, data, round)

    return sids, data


def read_response_data(dirstr):

    # Get data files
    dir   = './data/' + dirstr + '/'
    path  = dir + 'responses/'
    files = sorted(os.listdir(path))
    files = [file for file in files if file.startswith('p') and file.endswith('.csv')]

    # Manually marked as bad list
    bad_list = list(pd.read_csv(dir + '/remove_list.csv', header=1))

    # Read all the data
    sids, data = [], []
    for file in files:

        # Check if bad
        if file[1:5] in bad_list: continue

        # Subject id from filename
        sids.append(file[1:5])

        # Data from file, append sid
        df = pd.read_csv(path + file)
        df['sid'] = sids[-1]

        # The validation set (weirdly) has an extra NAN row.
        df = df.dropna(subset=['Question'])

        # Insert basic check that all subjects have same # of questions
        print('Subject '+ sids[-1] + ' has ' + str(df.shape[0]) + ' questions.')

        # Save to list
        data.append(df)

    # Merge data into single frame
    data = pd.concat(data).reset_index(drop = True)

    # Make sure we filtered properly
    assert not any([sid in bad_list for sid in sids])

    # Convert question number to an actual number
    data.quest_num = data.quest_num.apply(lambda x: int(x.split('/')[0]))

    return sids, data


#
# Attention check question 22 --> 13, 71 --> 42, 115 --> 71


def quality_check(sids, data, round):

    # Three check questions
    cqn = [22,71,115] if round == 1 else [13,42,71]
    check_1 = list(data.loc[data.quest_num == cqn[0],:].answer_num == 3)
    check_2 = list(data.loc[data.quest_num == cqn[1],:].answer_num == 1)
    check_3 = list(data.loc[data.quest_num == cqn[2],:].answer_num == 1)

    # Get run lengths for each subject
    runlens, check_4 = [], []
    for i, s in enumerate(sids):

        # Detect changes in answer
        df = pd.DataFrame()
        df['shifted'] = data[data.sid == s]['answer_num'].shift(1) != data[data.sid == s]['answer_num']

        # Cumulative sum of bools tells us which chunk (run) each answer falls in
        df['chunk'] = df['shifted'].cumsum()

        # Group them by run and count how many are in each run
        runlens.append( df.groupby('chunk').size().tolist() )

        # Check if any runs are longer than 10
        check_4.append( any([l < 10 for l in runlens[i]]) )

    # TODO: Add variance check back in?

    # List of pass/fail for each subject
    pass_check = []
    for i in range(len(check_1)):
        pass_check.append(check_1[i] and check_2[i] and check_3[i] and check_4[i])

    # Failed subject list
    failed = [sids[i] for i, val in enumerate(pass_check) if not val]

    # Display who failed
    if len(failed) == 0:
        print('No subjects failed attention checks.')
    else:
        print('Subjects failing checks:')
        print(failed)

    # Remove subjects failing checks from data
    for sid in failed:
        inds = data.index[data.sid == sid].tolist()
        data = data.drop(inds)
        data = data.reset_index(drop = True)

    # Remove subjects failing checks from sids list
    sids = [sid for sid in sids if sid not in failed]

    return sids, data


def get_qdata(data, sort = False):
    # Initialize the questions DataFrame with unique quest_num as the index.
    qnums  = np.unique(data.quest_num)
    nsubj  = len(np.unique(data.sid))

    qstats = pd.DataFrame(index = qnums, columns = ['mean', 'median', 'std', 'bin_use'])
    ansbyq = np.zeros([len(qnums), nsubj])
    qcdfs  = np.zeros([len(qnums), 4])
    qpdfs  = np.zeros([len(qnums), 4])

    # Calculate statistics for each unique question number
    for i in qnums:
        # Answers for this question
        answers = data.loc[data.quest_num == i, 'answer_num']
        
        # Aggregate answers into matrix
        ansbyq[i-1, :] = answers

        # Assign statistics
        qstats.loc[i, 'mean'   ] = np.mean(  answers)
        qstats.loc[i, 'median' ] = np.median(answers)
        qstats.loc[i, 'std'    ] = np.std(   answers)
        qstats.loc[i, 'var'    ] = np.var(   answers)

        # Number of bins used by participants
        qstats.loc[i, 'bin_use'] = len(np.unique(answers))

        # CDFs and PDFs
        qcdfs[i-1,:] = [np.sum(answers <= i+1)/nsubj for i in range(4)]
        qpdfs[i-1,:] = [np.sum(answers == i+1)/nsubj for i in range(4)]
        

    # Get indices of most to least informative
    inds = np.array(qstats['std'].sort_values(ascending = False).index) if sort else np.arange(len(qnums))+1
    
    # Sort everything by question standard deviation (loose informativeness)
    ansbyq = ansbyq[inds-1,:]
    qstats = qstats.loc[inds,:]
    qpdfs  = qpdfs[inds-1,:]
    qcdfs  = qpdfs[inds-1,:]

    # Questions in a list, in this order (loc is 0 indexed)
    qs = [Q[3:] for Q in data.loc[inds-1, 'Question']]

    return qs, ansbyq, qstats, qcdfs, qpdfs, inds


def get_pcs(ansbyq):
    # Covariance matrix

    # Get principal components
    d,v = np.linalg.eig(np.cov(ansbyq))

    # Convert to fraction variance explained
    d = np.real(d)
    varsum = np.sum(d)
    d = d/varsum

    d    = pd.DataFrame(d, columns = ['evals'])
    dinds = d['evals'].sort_values(ascending = False).index
    d = d.loc[dinds]
    v = v[:,dinds]

    return d, v



def run_sparse_PC_comparison(ansbyq, sample_frac = 1, resample = False, n_comp = 40, n_qsrm = 10, sparse = False):

    # Run baseline PCA
    pca = SparsePCA(n_components = n_comp) if sparse else PCA(n_components=n_comp)
    pca.fit(ansbyq.T)

    # Resample subject data
    if resample:
        # Number of subjects for data resampling
        n_subj = ansbyq.shape[1]
        n_subj = int(np.round(sample_frac*n_subj))

        # Resample the data
        inds = np.random.randint(72, size=[n_subj])
        respdata = ansbyq[:, inds]
    else:
        # Keep old
        respdata = ansbyq


    # Subset components from full PCA
    reduced_comps = pca.components_[:n_comp,:(150 - n_qsrm)]

    # Get PCA on pruned question set
    reduced_data     = ansbyq[0:(150 - n_qsrm),:]
    reduced_data_pca = sk.decomposition.SparsePCA(n_components = n_comp)
    reduced_data_pca.fit(reduced_data.T)

    # Components from new reduced-data PCA
    reduced_data_comps = reduced_data_pca.components_[:n_comp,:]

    # All components
    comps = np.concatenate([reduced_comps, reduced_data_comps], axis=0)

    # Correlation matrix
    plt.matshow(np.corrcoef(comps))



def sort_spca_components(spca, ansbyq, n_comp):

    # Question covariance
    qcov = np.cov(ansbyq)

    # Rank in order of variance explained
    ve = np.zeros(n_comp)
    for i in range(0, n_comp):
        ve[i] = np.matmul(spca.components_[i,:], np.matmul(qcov, spca.components_[i,:]))

    ve   = pd.DataFrame(ve, columns = ['ve'])
    inds = np.array(ve['ve'].sort_values(ascending = False).index)

    ve = np.array(ve).flatten()
    ve = ve[inds]

    spca.components_ = spca.components_[inds,:]

    return ve, spca


def run_pca_comparison(ansbyq):
    # Get baseline PCA results
    pca = PCA()
    pca.fit(ansbyq.T)

    # PCs we would like to keep
    pcs_to_keep = [1,2,5,6,8,9,15,26,37,39]

    # PCs we will try to keep
    n_top_pcs = 10
    high_load_inds = np.zeros([n_top_pcs,150])
    for i in range(0,n_top_pcs):
        high_load_inds[i,:] = abs(pca.components_[i]) > np.percentile(abs(pca.components_[i]),90)

    keep = np.sum(high_load_inds, axis = 0) > 0
    subsampled = ansbyq[keep,:]

    # New PCA results
    pcaB = PCA()
    pcaB.fit(subsampled.T)

    # Comparision
    ncomp = 15
    #corrs = np.corrcoef(np.concatenate([pca.components_[0:ncomp,keep], pcaB.components_[0:ncomp,:]]))
    allvecs = np.concatenate([pca.components_[0:ncomp,keep], pcaB.components_[0:ncomp,:]])
    ips = np.matmul(allvecs, allvecs.T)
    plt.matshow(ips,aspect='auto')

    return keep

    # For Noham
    # n_pcs_to_keep = 11
    # high_load_inds = np.zeros([n_pcs_to_keep,150])
    # for i, n in enumerate(pcs_to_keep):
    #     high_load_inds[i,:] = abs(pca.components_[n]) > np.percentile(abs(pca.components_[n]),95)

    # keep = np.sum(high_load_inds, axis = 0) > 0

    # subsampled = ansbyq[keep,:]

    # # New PCA results
    # pcaB = PCA()
    # pcaB.fit(subsampled.T)

    # # Comparision
    # corrs = np.corrcoef(np.concatenate([pca.components_[:,keep], -pcaB.components_]))
    # plt.matshow(corrs[pcs_to_keep,:],aspect='auto')


