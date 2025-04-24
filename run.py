from configs import *

# Read subject data
sids, data = import_data(subj_data_dir)

# Chunk question data, get basic stats
qs, ansbyq, qstats, qcdfs, qpdfs, inds = get_qdata(sids, data, sort=False)

# Manually generated list of reverse coded questions
reverse_coded = np.array([
    1 ,3 , 5 ,6 , 11,12, 14,15, 18,20,
    22,23, 27,28, 30,31, 35,36, 38,40,
    41,44, 46,48, 51,52, 54,55, 57,59,
    62,64, 67,68, 71,72, 75,76, 78,80,
    83,84, 86,87, 90,92 ,94,96, 98,100]) -1

# Flip the reverse coded questions
recoded = flip_reverse_coded(ansbyq, reverse_coded)

# Plot the data
plot_means(qstats)
plot_stds(qstats)

# Generate a copy with NaNs removed
nanlocs = np.any(np.isnan(recoded),axis=1)
reduced = recoded[~nanlocs,:]

# Plot correlation matrix
plot_corrs(reduced)

# PCA of questionnaire data
pca = PCA()
pca.fit(reduced.T)

# Plot PCA results 
plot_pca_variance_explained( pca)
plot_pca_cumulative_variance(pca)
plot_pca_components(pca)

# Get abbreviated questions for plots
short_qs = [' '.join(Q.split(' ')[0:10]) for Q in qs]

# Plot top and bottom 40 most informative PDFs
plot_response_dists(qpdfs, np.arange(0,40)  , 'Response PDFs (Top 40)', short_qs)
plot_response_dists(qpdfs, np.arange(-41,-1), 'Response PDFs (Bottom 40)', short_qs)

# Plot variances in original grid (will want to make these "collapsible" across domain & construct too)
std_heatmap, qcoords = create_construct_domain_heatmap(qstats['std'].values)
plot_construct_domain_heatmap(std_heatmap)