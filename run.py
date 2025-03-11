from configs import *

def response_data_analysis(qs, ansbyq, qstats, qpdfs):
    # Plot answer means and standard deviations
    plot_means(qstats)
    plot_stds(qstats)

    # PCA of questionnaire data
    pca = PCA()
    pca.fit(ansbyq.T)

    # Plot PCA results 
    plot_pca_variance_explained( pca)
    plot_pca_cumulative_variance(pca)

    # Get abbreviated questions for plots
    short_qs = [' '.join(Q.split(' ')[0:10]) for Q in qs]

    # Plot top and bottom 40 most informative PDFs
    plot_response_dists(qpdfs, np.arange(0,40)   , 'Response PDFs (Top 40)', short_qs)
    plot_response_dists(qpdfs, np.arange(-41,-1), 'Response PDFs (Bottom 40)', short_qs)


# Get subject response data
sids_1, data_1 = import_data('round-1', round = 1)

# Chunk question data, get basic stats
qs_1, ansbyq_1, qstats_1, qcdfs_1, qpdfs_1, inds_1 = get_qdata(data_1)

# Analyse the data
response_data_analysis(qs_1, ansbyq_1, qstats_1, qpdfs_1)

# Compare baseline and potential reduced questionnaire PCs
keep = run_pca_subsample_analysis(ansbyq_1)

# Do all the same things for the second round of response data
sids_2, data_2 = import_data('round-2', round = 2)
qs_2, ansbyq_2, qstats_2, qcdfs_2, qpdfs_2, inds_2 = get_qdata(data_2)
response_data_analysis(qs_2, ansbyq_2, qstats_2, qpdfs_2)

# Figure out which questions in round 1 matched round 2. Zero indexed!
match_pairs = []
for i in range(94): 
    q2  = data_2.loc[data_2['sid'] == sids_2[0]]['Question'][i]
    q1s = data_1.loc[data_1['sid'] == sids_1[0]]['Question']
    match = np.where( q2 == q1s)[0][0]

    match_pairs.append([i, match])

    #print(' ')
    #print(q2)
    #print(q1s[match])
    #print(str(i) + ' ' + str(match))
    #input("Press Enter to continue...")

match_pairs = np.array(match_pairs)
match_inds  = match_pairs[:,1]

assert [qs_1[i] for i in match_inds] == qs_2

# Concatenate all subject data for shared questions for comparison
ansbyq_3 = np.concatenate([ansbyq_1[match_inds,:],ansbyq_2],axis=1)

# Compare PCs from first, subsampled, second round, and combined data
ncs = 10

# Original PCs
pca_1 = PCA(n_components = ncs)
pca_1.fit(ansbyq_1.T)

# PCs of subsampled data
pca_2 = PCA(n_components = ncs)
pca_2.fit(ansbyq_1[match_inds,:].T)

# PCs of round 2
pca_3 = PCA(n_components = ncs)
pca_3.fit(ansbyq_2.T)

# PCs on all data
pca_4 = PCA(n_components = ncs); 
pca_4.fit(ansbyq_3.T)

# Plot inner products of round 1 matched subset PCs and round 2 PCs
cinds = np.arange(0,ncs)
stacked = [pca_2.components_[cinds,:], pca_3.components_[cinds,:]]
stacked = np.concatenate(stacked)
ips = np.matmul(stacked, stacked.T)
plt.matshow(ips, aspect = 'auto')

# Plot inner products of truncated round 1 PCs and round 1 matched subset PCs
stacked = np.concatenate([pca_1.components_[cinds,:][:,match_inds], pca_2.components_[cinds,:]])
ips = np.matmul(stacked, stacked.T)
plt.matshow(ips,aspect='auto')

# Plot inner products of truncated round 1 PCs and round 2 PCs
stacked = np.concatenate([pca_1.components_[cinds,:][:,match_inds], -pca_3.components_[cinds,:]])
ips = np.matmul(stacked, stacked.T)
plt.matshow(ips,aspect='auto')

# Plot inner products of truncated, subsetted, and round 2, and joint
stacked = np.concatenate([pca_1.components_[cinds,:][:,match_inds], pca_2.components_[cinds,:], pca_3.components_[cinds,:], pca_4.components_[cinds,:]])
ips = np.matmul(stacked, stacked.T)
plt.matshow(ips,aspect='auto')

