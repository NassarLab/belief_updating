import numpy as np

#def factor_analysis_subsets(ansbyq, qs, qstats):
plt.close('all')
use_fa = True

# random exclusions
exclusions = np.random.choice(range(0,76), size=10, replace=False)
keep = np.setdiff1d(range(0,76), exclusions)

nqs = 6
#factor_data = pd.read_excel('/home/dan//Downloads/sorted_by_var_by_loading_paf_oblimin.xlsx')

newbyq = ansbyq.copy()
newbyq = newbyq[:, keep]

if use_fa:
    fa, loadings = factor_analysis(newbyq, qs=qs, qvar=qstats['var'], rotation='oblimin', n_factors=3, method='principal', sklearn=False)
else:
    pca = PCA()
    pca.fit(newbyq.T)
    loadings = pd.DataFrame(pca.components_[0:3].T, columns=['factor' + str(i+1) for i in range(3)])
    loadings['variance'] = qstats['var']
    loadings['item'] = qs

factor_data = loadings.sort_values(by='variance', ascending=False)
factor_data.reset_index(inplace=True, drop=True)

lower_std = np.percentile(qstats['var'], 25)
upper_std = np.percentile(qstats['var'], 100)

nrows = len(factor_data)
ratios = np.zeros(nrows)
f1_ratios = np.zeros(nrows)
f2_ratios = np.zeros(nrows)
f3_ratios = np.zeros(nrows)
for i in range(0, nrows):
    row = factor_data.iloc[i].values[0:3].astype(float)
    max = np.max(np.abs(row))
    sum = np.sum(np.abs(row))
    ratios[i] = max/sum

    f1_ratios[i] = np.abs(row)[0] / np.sum(np.abs(row))
    f2_ratios[i] = np.abs(row)[1] / np.sum(np.abs(row))
    f3_ratios[i] = np.abs(row)[2] / np.sum(np.abs(row))

sort_inds_ratios = np.flip(np.argsort(ratios))

sort_inds_f1_ratios = np.flip(np.argsort(f1_ratios))[0:30]
sort_inds_f2_ratios = np.flip(np.argsort(f2_ratios))[0:30]
sort_inds_f3_ratios = np.flip(np.argsort(f3_ratios))[0:30]

#print('Top 10 item ratios for each factor:')
#print(f1_ratios[sort_inds_f1_ratios][0:10])
#print(f2_ratios[sort_inds_f2_ratios][0:10])
#print(f3_ratios[sort_inds_f3_ratios][0:10])


f1_inds = np.array([i for i in sort_inds_f1_ratios if (factor_data['variance'][i] > lower_std and factor_data['variance'][i] < upper_std)])
f2_inds = np.array([i for i in sort_inds_f2_ratios if (factor_data['variance'][i] > lower_std and factor_data['variance'][i] < upper_std)])
f3_inds = np.array([i for i in sort_inds_f3_ratios if (factor_data['variance'][i] > lower_std and factor_data['variance'][i] < upper_std)])


# print(factor_data.iloc[f1_inds][0:nqs])
# print('')
# print(factor_data.iloc[f2_inds][0:nqs])
# print('')
# print(factor_data.iloc[f3_inds][0:nqs])

###
keep = np.concatenate([f1_inds[0:nqs], f2_inds[0:nqs], f3_inds[0:nqs]])

f1_indices_orig = [i for i,q in enumerate(qs) if q in factor_data.iloc[f1_inds]['item'].tolist()]
f2_indices_orig = [i for i,q in enumerate(qs) if q in factor_data.iloc[f2_inds]['item'].tolist()]
f3_indices_orig = [i for i,q in enumerate(qs) if q in factor_data.iloc[f3_inds]['item'].tolist()]

keep_orig = np.concatenate([f1_indices_orig[0:nqs], f2_indices_orig[0:nqs], f3_indices_orig[0:nqs]])

pca = PCA()
pca.fit(zbyq[keep_orig,:].T)


def plot_everything():
    # Plot PCA results 
    plot_pca_variance_explained(pca, title='Item Level PCA Variance Explained')
    plot_pca_components(pca, xlabel='Item', title='Item Level PCA Components', comps=[0,1,2])


    # Get big correlation matrix
    qqcorrs = np.corrcoef(ansbyq)

    # Plot correlation matrix
    plot_qq_corrs(ansbyq, keep_orig)


    plot_corrs_by_means(ansbyq, forward_coded, reverse_coded, keep_orig)
    plot_stds_by_means(ansbyq, forward_coded, reverse_coded, keep_orig)
    plot_corrs_by_stds(ansbyq, forward_coded, reverse_coded, keep_orig)

    plot_qstds(qstats, forward_coded, reverse_coded, keep_orig)


    for q in [qs[i] for i in f1_indices_orig[0:6]]:print(q)
    print('')
    for q in [qs[i] for i in f2_indices_orig[0:6]]:print(q)
    print('')
    for q in [qs[i] for i in f3_indices_orig[0:6]]:print(q)