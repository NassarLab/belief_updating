import numpy   as np
import pandas  as pd
import sklearn as sk
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA, SparsePCA, FactorAnalysis
from factor_analyzer import FactorAnalyzer


def get_qdata(sids, data, sort = False, trim=True, flip=True):
    """
    Takes the stacked subject data and extracts question information.

    If sort is true, the questions are sorted by standard deviation.
    If trim is true, attention checks are removed.
    If flip is true, reverse coded questions are flipped.

    Returns:
    qs:      List of questions
    ansbyq:  Answer matrix (questions x subjects)
    qstats:  Question statistics (mean, median, std, var, bin_use)
    qcdfs:   Question CDFs (questions x Likert options)
    qpdfs:   Question PDFs (questions x Likert options)
    inds:    Indices of questions sorted from high to low variance
    """
    print('Getting question x answer data and related stats ...')

    # Get question numbers, potentially remove attention checks
    qnums = np.unique(data['original_index'])
    if trim: qnums = qnums[~np.isin(qnums, [101, 102, 103, 104])]

    # Number of subjects, questions 
    nsubj = len(np.unique(data.sid))
    nqs   = len(qnums)

    # Indices of reverse coded questions
    reverse_coded = np.array([
        1 ,3 , 5 ,6 , 11,12, 14,15, 18,20,
        22,23, 27,28, 30,31, 35,36, 38,40,
        41,44, 46,48, 51,52, 54,55, 57,59,
        62,64, 67,68, 71,72, 75,76, 78,80,
        83,84, 86,87, 90,92 ,94,96, 98,100]) -1

    # Question statistics, answer matrix, CDFs and PDFs
    ansbyq = np.zeros([nqs, nsubj])
    qstats = np.zeros([nqs, 5])
    qcdfs  = np.zeros([nqs, 5])
    qpdfs  = np.zeros([nqs, 5])

    # Initialize list of questions to reorder if sort is true.
    qs = []

    # Calculate statistics for each unique question number
    for qn in qnums.astype(int):
        # Questions are 1-indexed in data, should be 0-indexed in arrays here
        i = qn - 1
        
        # Answers for this question
        q_matches = data.loc[data['original_index'] == qn, ['answer_num','sid']]
        q_answers = np.zeros(nsubj)*np.nan
        
        # Check whether there are too many answers for some reason
        if len(q_matches) != nsubj:
            print('Warning: Question ' + str(qn) + ' has ' + str(len(q_matches)) + ' answers, not ' + str(nsubj) + '.')

        # Get each subject's set of answers to this question
        for j, sid in enumerate(sids):
            sq_answers = q_matches.loc[q_matches.sid == sid, 'answer_num']

            # If there aren't any, set to NaN
            if len(sq_answers) == 0:
                print('Warning: Subject ' + sid + ' has no answer for question ' + str(qn) + '.')
                q_answers[j] = np.nan

            # If there is only one (as there should be), use it
            if len(sq_answers) == 1:
                q_answers[j] = sq_answers.iloc[0]

            # If there is more than one, take the last
            if len(sq_answers) > 1:
                print('Warning: Subject ' + sid + ' has multiple answers for question ' + str(qn) + '.')
                q_answers[j] = sq_answers.iloc[-1]

        # Check if we need to flip the answers
        if flip and (i in reverse_coded): q_answers = 6 - q_answers            

        # Save row to return array
        ansbyq[i, :] = q_answers

        # Get statistics
        qstats[i,:] = [np.mean(q_answers), np.median(q_answers), np.std(q_answers), np.var(q_answers), len(np.unique(q_answers))]

        # CDFs and PDFs
        qcdfs[i,:] = [np.sum(q_answers <= j)/nsubj for j in range(1,6)]
        qpdfs[i,:] = [np.sum(q_answers == j)/nsubj for j in range(1,6)]
        
        qs.append(data.loc[data['original_index'] == qn, ['Question']].values[0][0])

    # Convert qstats to DataFrame
    qstats = pd.DataFrame(qstats, columns = ['mean', 'median', 'std', 'var', 'bin_use'])
    qstats['question'] = qs
    qstats['qnum']  = qnums

    # Sort indices of most to least informative if requested
    inds = np.flip(np.argsort(qstats['std'].astype(float).values))
    
    # Sort everything by question standard deviation if requested
    if sort:
        qs = [qs[i] for i in inds]
        ansbyq = ansbyq[inds,:]
        qstats = qstats.loc[inds,:].reset_index(drop = True)
        qpdfs  = qpdfs[inds,:]
        qcdfs  = qpdfs[inds,:]

    return qs, ansbyq, qstats, qcdfs, qpdfs, inds



def run_sparse_PC_comparison(ansbyq, sample_frac = 1, resample = False, n_comp = 40, n_qsrm = 10, sparse = False):
    """
    Old function to compare PCs from multiple rounds of data collection, needs updating.
    """

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
    """
    Old function used in compareing PCs from multiple rounds of data collection, needs updating.
    """

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


def run_pca_subsample_analysis(ansbyq, ):
    """
    Old function used in comparing subsampled PCs and PCs from subsampled data, needs updating.
    """
    # Get baseline PCA results
    pca = PCA()
    pca.fit(ansbyq.T)

    # PCs we would like to keep
    # pcs_to_keep = [1,2,5,6,8,9,15,26,37,39]

    # Try to keep top PCs using most loaded questions
    n_top_pcs = 10
    high_load_inds = np.zeros([n_top_pcs,150])
    for i in range(0,n_top_pcs):
        # Check which questions have high loadings
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


def factor_analysis(zbyq, qs, qvar=None, rotation='oblimin', n_factors=3, method='principal', sklearn=False):
    
    df = pd.DataFrame(zbyq.T, columns=['Q' + str(i+1) for i in range(100)])

    # Perform factor analysis
    if sklearn:
        fa = FactorAnalysis(n_components=n_factors, rotation=rotation)
        fa.fit(df)
        loadings = pd.DataFrame(fa.components_.T, columns=['factor' + str(i+1) for i in range(n_factors)])
        print('Sklearn factor analysis does not support different methods.')
    else:
        fa = FactorAnalyzer(n_factors=n_factors, rotation=rotation, method=method)
        fa.fit(df)
        loadings = pd.DataFrame(fa.loadings_, columns=['factor' + str(i+1) for i in range(n_factors)])

    # Save any additional fields
    if qvar is not None:
        loadings['variance'] = qvar
    loadings['item'] = qs

    return fa, loadings

# Compute running variance in sliding scale over 10 questions for each participant
def get_running_average_response_stats(sids, data):
    window = np.ones(10)/10
    running_avg = np.full((len(sids), len(qs)-10+1), np.nan)
    running_std = np.full((len(sids), len(qs)-10+1), np.nan)
    for i in range(len(sids)):
        answers = data.loc[data['sid']==sids[i]]['answer_num'].values[0:100]
        num_qs = len(answers)
        if num_qs < 91:
            print('Participant {} has {} questions'.format(sids[i], num_qs))
            continue
        try:
            running_avg[i,:] = np.convolve(answers, window, mode='valid')
            for j in range(len(qs)-10+1):
                running_std[i,j] = np.std(answers[j:(j+10)])
        except:
            continue

    return running_avg, running_std