
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