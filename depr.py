
# Get cut-down PCA results
pca50 = PCA()
pca50.fit(ansbyq[0:50,:].T)

# Get cut-down PCA results
pca100 = PCA()
pca100.fit(ansbyq[0:100,:].T)


# Get baseline PCA results
spca = SparsePCA(alpha = 1, n_components=40, verbose = 2, max_iter=1000)
spca.fit(ansbyq.T)

# Get cut-down PCA results
spca50 = SparsePCA(alpha = 1, n_components=40, verbose = 2, max_iter=1000)
spca50.fit(ansbyq[0:50,:].T)

# Get cut-down PCA results
spca100 = SparsePCA(alpha = 1, n_components=40, verbose = 2, max_iter=1000)
spca100.fit(ansbyq[0:100,:].T)

ve   , spca    = sort_spca_components(spca   , ansbyq         , 40)
ve100, spca100 = sort_spca_components(spca100, ansbyq[0:100,:], 40)
ve50 , spca50  = sort_spca_components(spca50 , ansbyq[0:50 ,:], 40)

plot_response_dists(spca, spca50, spca100, ncomp = 40)




from sklearn.decomposition import FastICA as ICA

ncs = 10
pca_1 = decomp(n_components = ncs)
pca_1.fit(ansbyq_1[match_inds,:].T)

pca_2 = decomp(n_components = ncs)
pca_2.fit(ansbyq_2.T)

cinds = np.arange(0,ncs)
allvecs = [pca_1.components_[cinds,:], pca_2.components_[cinds,:]]
allvecs = np.concatenate(allvecs)
#plt.matshow(corrs[0:10,:],aspect='auto')
plt.matshow(corrs,aspect='auto')





#for i in range(0,10):
plt.figure()
plt.plot(  pca_2.components_[i,:])
plt.plot(- pca_3.components_[i,:])
#plt.plot(  pca_4.components_[i,:])
plt.title(f"PC {i:1d}")


dem_2 = pd.read_csv('./data/round-2/demographic_information.csv')

keep = []
for i in range(0,94):
    print(' ')
    print(qs_2[i])
    inp = input("Press Enter to continue...")
    if inp == 'y': keep.append(qs_2[i]) 



plt.figure()
vals, cdf = get_cdf(dem_1['Age'])
plt.plot(vals, cdf)
vals, cdf = get_cdf(dem_2['Age'])
plt.plot(vals, cdf)

plt.figure()
vals, cdf = get_cdf(dem_1['Time taken'])
plt.plot(vals, cdf)
vals, cdf = get_cdf(dem_2['Time taken'])
plt.plot(vals, cdf)
vals, cdf = get_cdf(dem_1['Time taken']*(94/150))
plt.plot(vals, cdf)




plt.figure()
vals, cdf = get_cdf(dem_1['Sex'])
plt.plot(vals, cdf)
vals, cdf = get_cdf(dem_2['Sex'])
plt.plot(vals, cdf)


plt.figure()
vals, cdf = get_cdf(dem_1['Total approvals'])
plt.plot(vals, cdf)
vals, cdf = get_cdf(dem_2['Total approvals'])
plt.plot(vals, cdf)


plt.figure()
vals, cdf = get_cdf(dem_1['Ethnicity simplified'])
plt.plot(vals, cdf)
vals, cdf = get_cdf(dem_2['Ethnicity simplified'])
plt.plot(vals, cdf)


plt.figure();
plt.plot(dem_1.loc[dem_1['Sex'] == 'Female']['Age'],'o')
plt.plot(dem_1.loc[dem_1['Sex'] == 'Male']['Age'],'o')

plt.figure();
plt.plot(dem_2.loc[dem_2['Sex'] == 'Female']['Age'],'o')
plt.plot(dem_2.loc[dem_2['Sex'] == 'Male']['Age'],'o')



plt.figure(); 
plt.plot(np.arange(94), qstats_1.loc[match_pairs[:,1]+1,:]['mean'],'o')
plt.plot(np.arange(94), qstats_2['mean'],'o')



plt.figure(); 
plt.plot(np.arange(94), qstats_1.loc[match_pairs[:,1]+1,:]['std'],'o')
plt.plot(np.arange(94), qstats_2['std'],'o')


def previous_rounds_runscript():
    # Get subject response data
    sids_1, data_1 = import_data('take-2', round = 1)

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



def get_cdf(data, type = 'cdf'):
    vals = np.unique(data)
    cnts = np.array([sum(item == data) for item in vals])
    csum = np.cumsum(cnts)
    cdf  = csum/csum[-1]
    pdf  = cnts/sum(cnts)

    df = cdf if type == 'cdf' else pdf

    return vals, df




# def plot_response_dists(pca, pca50, pca100, ncomp):
    
#     # Get correlations
#     corrs = np.corrcoef(np.concatenate([pca.components_[0:ncomp,0:50], pca100.components_[0:ncomp,0:50], pca50.components_[0:ncomp,0:50] ]))

#     # Plot response DFS
#     fig, ax = plt.subplots(figsize = [6,6])
#     cax = ax.matshow(corrs, aspect='auto')

#     plt.title('PC Correlations')

#     #ax.set_xticks([i for i in range(0,15)])
#     #ax.set_xticklabels(['C' + str(i) + str(j) for i in range(0,3) for j in range(0,5)])

#     #ax.set_yticks([i for i in range(0,15)])
#     #ax.set_yticklabels(['C' + str(i) + str(j) for i in range(0,3) for j in range(0,5)])

#     plt.tight_layout()

#     return corrs
    