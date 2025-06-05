from configs import *

# Read subject data
sids, data = read_response_data(subj_data_dir)

# Chunk question data, get basic stats
qs, ansbyq, qstats, qcdfs, qpdfs, inds = get_qdata(sids, data, sort=False, trim=True, flip=True)

# Get construct, domain, coding, and grid coordinates by question
qcoords = get_general_question_coordinates()

# Get aggregated scores by construct, domain, and reverse coding
zbyq, agg_cstr, agg_cstr_dom, agg_cstr_dom_rev = ans_to_scores(ansbyq, qcoords)

# Get abbreviated questions for plots
short_qs = [' '.join(Q.split(' ')[0:5]) for Q in qs]

# Reverse coded questions
reverse_coded = np.array([
    1 ,3 , 5 ,6 , 11,12, 14,15, 18,20,
    22,23, 27,28, 30,31, 35,36, 38,40,
    41,44, 46,48, 51,52, 54,55, 57,59,
    62,64, 67,68, 71,72, 75,76, 78,80,
    83,84, 86,87, 90,92 ,94,96, 98,100]) -1

forward_coded = np.setdiff1d(np.arange(100), reverse_coded)


# ---- Plot means, variances, and bin use in original grid ---- #
mean_heatmap = create_construct_domain_item_heatmap(qstats['mean'].values, qcoords)
plot_construct_domain_item_heatmap(mean_heatmap, vmin=1, vmax=5, title='Mean Values')

unif_std = 1.414
std_heatmap = create_construct_domain_item_heatmap(qstats['std'].values/unif_std, qcoords)
plot_construct_domain_item_heatmap(std_heatmap, vmin=0, vmax=1, title='Normalized Standard Deviations')

bin_heatmap = create_construct_domain_item_heatmap(qstats['bin_use'].values, qcoords)
plot_construct_domain_item_heatmap(bin_heatmap, vmin=1, vmax=5, title='Bin Use')

low_std  = qstats['std'] <= np.percentile(qstats['std'], 25)
low_hmap = create_construct_domain_item_heatmap(low_std.values, qcoords)
plot_construct_domain_item_heatmap(low_hmap, title='Low Standard Deviation Items')

# ---- Plot item CDFS --- #
plot_response_dists(qcdfs, lb=0 , ub=33 , title='Response CDFs 1-33'  , qs=short_qs)
plot_response_dists(qcdfs, lb=33, ub=66 , title='Response PDFs 34-66' , qs=short_qs)
plot_response_dists(qcdfs, lb=67, ub=100, title='Response PDFs 67-100', qs=short_qs)


# ---- Plot raw item correlation matrix ---- #
plot_response_corrs_anotated(ansbyq)


# ---- Plot PCA of the data ---- #
pca = PCA()
pca.fit(zbyq.T)
ncomps = 3
plot_pca_variance_explained(pca, title='Item Level PCA Variance Explained')
plot_pca_components(pca, xlabel='Item', title='Item Level PCA Components', comps=range(0,ncomps))
plot_component_heatmaps(pca, qcoords, ncomps)

# ---- Plot PCA of construct scores ---- #
pca = PCA()
pca.fit(agg_cstr)
constructs = ['Det', 'Inc', 'Alt', 'Inf', 'Opn']
plot_pca_components(pca, xlabel='Construct', title='Construct Score PCs', comps=[0,1,3], xticklabels=constructs)
plot_pca_variance_explained(pca, 'Between-Construct VEs')

# ---- Construct score SDs ---- #
plot_construct_score_stds(agg_cstr, constructs)

# ------- For each construct, perform a PCA on the scores to get a factor dominance ----- #
plot_construct_PCAs(constructs, qcoords, zbyq)

# ---- Plot Aggregated Correlations ---- #
plot_constr_score_corrs(agg_cstr)
plot_domain_score_corrs(agg_cstr_dom)
plot_reverse_coding_score_corrs(agg_cstr_dom_rev)
plot_domain_score_corrs_separated(agg_cstr_dom)

# ---- Plot Domain Aggregated Statistics ---- #
cstr_dom_stds = np.std(agg_cstr_dom, axis=0)
plot_construct_domain_heatmap(cstr_dom_stds, qcoords, title='Construct-Domain Standard Deviations')


# ---- Plot PCA of domain scores ---- #
plot_domain_score_pca(agg_cstr_dom, qcoords)


# -------------------- Item removal, low SDs and no incrementalism -------------------- #
plot_qstds(qstats, forward_coded, reverse_coded)

# Drop low STD items and those from incrementalism construct
low_std  = qstats['std'] >= np.percentile(qstats['std'], 50)
drop_low_idx = np.where(low_std)[0]
drop_inc_idx = np.where(qcoords['cnum'] == 1)[0]
drop_idx = np.concatenate([drop_low_idx, drop_inc_idx])

plot_pca_of_subsetted_questions(ansbyq, qcoords, drop_idx, ncomps=3)

# ----------- Plot participant total stds --------------#
plot_participant_grand_stds(ansbyq)


# ------------- Plot running averages and standard deviations -------------- #

# Compute running variance in sliding scale over 10 questions for each participant
running_avg, running_std = get_running_average_response_stats(sids, data)

# Plot running averages and their trends over survey time
plot_running_average_responses(running_avg)
plot_running_average_trends(running_avg, sids)

# Plot running standard deviations and their trends over survey time
plot_running_std(running_std)
plot_running_std_trends(running_std, sids)


# -------------- Plot question correlations, means, and standard deviations relative to one another
plot_corrs_by_means(ansbyq, forward_coded, reverse_coded)
plot_stds_by_means(ansbyq, forward_coded, reverse_coded)
plot_corrs_by_stds(ansbyq, forward_coded, reverse_coded)

# Plot question SD distribution (sorted)
plot_qstds(qstats, forward_coded, reverse_coded)


# -------------- Top factor analyses ------------------------#