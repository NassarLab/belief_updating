
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