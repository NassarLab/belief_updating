import numpy as np
import matplotlib.pyplot as plt

def plot_response_dists(dfs, lb, ub, title, qs):
    # Indices from bounds
    inds = np.arange(lb, ub)

    # Plot response DFS
    fig, ax = plt.subplots(figsize = [5.47, 7.51])
    cax = ax.matshow(dfs[inds,:], aspect='auto', vmin=0, vmax=1)
    plt.title(title)

    # Set questions as y-axis tick labels
    ax.yaxis.set_ticks_position('right')
    ax.set_yticks(inds-inds[0])
    ax.set_yticklabels([qs[i] for i in inds])

    # Set agreement as x-axis tick labels
    ax.set_xticks([0,1,2,3,4])
    ax.set_xticklabels(['Strongly Disagree','Slightly Disagree', 'Neutral', 'Slightly Agree', 'Strongly Agree'], rotation = 30, ha='left')

    plt.tight_layout()


# Plot standard deviations
def plot_stds(qstats):
    qnums = np.arange(0, qstats.shape[0])
    plt.figure(figsize = [4,4])
    plt.plot(qnums, qstats['std'].values, 'o')
    plt.title('Answer Standard Deviations')
    plt.xlabel('Sorted Question Number')
    plt.ylabel('Standard Deviation')
    plt.tight_layout()


# Plot means
def plot_means(qstats):
    qnums = np.arange(0, qstats.shape[0])
    plt.figure(figsize = [4,4])
    plt.plot(qnums, qstats['mean'].values, 'o')
    plt.title('Answer Means')
    plt.xlabel('Sorted Question Number')
    plt.ylabel('Mean')
    plt.tight_layout()

def plot_pca_variance_explained(pca, title='PCA Scree Plot'):
    # Plot relative variance explained
    plt.figure(figsize = [4,4])
    plt.plot(pca.explained_variance_ratio_, '-o')
    plt.title(title)
    plt.xlabel('Component')
    plt.ylabel('Fraction Variance Explained')
    #plt.grid('on')
    plt.tight_layout()

def plot_pca_cumulative_variance(pca, title='PCA Cumulative Variance Explained'):
    # Plot cumulative variance explained
    plt.figure(figsize = [4,4])
    plt.plot(np.cumsum(pca.explained_variance_ratio_),'-o')
    plt.title(title)
    plt.xlabel('Component')
    plt.ylabel('Cumulative Relative Variance')
    plt.grid('on')
    plt.tight_layout()

def plot_pca_components(pca, xlabel, title, comps=[0,1], xticklabels=None):
    # Plot cumulative variance explained
    plt.figure(figsize = [4,4])
    plt.plot(pca.components_[comps,:].T,'-o')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel('Loading')
    if xticklabels is not None:
        plt.xticks(np.arange(0, len(xticklabels)), xticklabels)
    #plt.grid('on')
    plt.axhline(y=0, color='k', linestyle='--')
    plt.tight_layout()

def plot_response_corrs_annotated(data):
    # Get correlations
    corrs = np.corrcoef(data)

    # Plot the data
    plt.figure(figsize=[6.12, 5.35])
    plt.imshow(corrs, aspect='auto', vmin=-1, vmax=1, interpolation='none')
    plt.gca().tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
    plt.title('Answer Correlation Matrix')
    plt.xlabel('Question')
    plt.ylabel('Question')
    plt.colorbar(label='Correlation Coefficient')

    # Insert lines every 20th row and column
    for i in range(1, 5):
        plt.axhline(y=i*20-0.5, color='k', linewidth=2)
        plt.axvline(x=i*20-0.5, color='k', linewidth=2)
    
    # Label groups
    plt.xticks(np.arange(0, 100, 20)-0.5+10, ['Det.', 'Incr.', 'Alt.', 'Info.', 'Open.'])
    plt.yticks(np.arange(0, 100, 20)-0.5+10, ['Det.', 'Incr.', 'Alt.', 'Info.', 'Open.'])

    plt.tight_layout()


def plot_construct_domain_item_heatmap(heatmap, vmin=None, vmax=None, title='Construct-Domain Heatmap'):
    """
    Plots a heatmap on the excel grid. Heatmap should be on this grid already.
    """
    # Plot the data
    plt.figure(figsize = [3.9, 6.7])
    plt.imshow(heatmap, aspect='auto', interpolation='None')
    plt.title(title)
    plt.gca().tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)

    if vmin is not None:
        plt.clim(vmin, vmax)

    # Insert horizontal lines at every 4th row, vertical at every column
    for i in range(1, 5):
        plt.axhline(y=i*4-0.5, color='k', linewidth=2)
        plt.axvline(x=i-0.5, color='k', linewidth=2)

    # Change x and y ticks
    plt.xticks(np.arange(5), ['General', 'Social', 'Self', 'Policy', 'News'])
    plt.yticks(np.arange(1.5, 20, 4), ['Detail', 'Incr.', 'Alt.', 'Info', 'Open'])

    # Insert value into each cell
    nrows, ncols = heatmap.shape
    for i in range(nrows):
        for j in range(ncols):
            plt.text(j, i, str(round(heatmap[i,j], 2)), ha='center', va='center')

    plt.tight_layout()

def plot_construct_domain_heatmap(heatmap, qcoords, vmin=None, vmax=None, title='Construct-Domain Heatmap'):
    """
    Plots a heatmap on the excel grid. Heatmap should be on this grid already.
    """
    # Regrid to item-level
    new_heatmap = np.zeros([5*4, 5])
    for qnum in range(100):
        i,j = qcoords.iloc[qnum][['i', 'j']]
        cnum, dnum = qcoords.iloc[qnum][['cnum', 'dnum']]

        # Regrid the data, save the q-coordinates
        new_heatmap[i,j] = heatmap[cnum, dnum]


    # Plot the data
    plt.figure(figsize = [3.9, 6.7])
    plt.imshow(new_heatmap, aspect='auto', interpolation='None')
    plt.title(title)
    plt.gca().tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)

    if vmin is not None:
        plt.clim(vmin, vmax)

    # Insert horizontal lines at every 4th row, vertical at every column
    for i in range(1, 5):
        plt.axhline(y=i*4-0.5, color='k', linewidth=2)
        plt.axvline(x=i-0.5, color='k', linewidth=2)

    # Change x and y ticks
    plt.xticks(np.arange(5), ['General', 'Social', 'Self', 'Policy', 'News'])
    plt.yticks(np.arange(1.5, 20, 4), ['Detail', 'Incr.', 'Alt.', 'Info', 'Open'])

    # Insert value into each cell
    nrows, ncols = new_heatmap.shape
    for i in range(nrows):
        for j in range(ncols):
            if i % 4 == 0:
                plt.text(j, i+2-0.5, str(round(new_heatmap[i,j], 2)), ha='center', va='center')

    plt.tight_layout()





def plot_constr_score_corrs(agg_cstr, colorbar=False):
    corrs = np.corrcoef(agg_cstr.T)
    plt.figure(figsize = [4.1, 4.2])
    plt.imshow(corrs, aspect='auto', vmin=-1, vmax=1, interpolation='none')
    plt.xticks(np.arange(5), ['Detail', 'Incr.', 'Alt.', 'Info', 'Open'])
    plt.yticks(np.arange(5), ['Detail', 'Incr.', 'Alt.', 'Info', 'Open'])
    plt.title('Construct Score Correlation Matrix')
    plt.xlabel('Construct')
    plt.ylabel('Construct')

    # Label each cell with numeric value
    nrows, ncols = corrs.shape
    for i in range(nrows):
        for j in range(ncols):
            plt.text(j, i, str(round(corrs[i,j], 2)), ha='center', va='center')

    if colorbar: plt.colorbar()
    plt.tight_layout()


def plot_domain_score_corrs(agg_cstr_dom):
    ns, nc, nd = agg_cstr_dom.shape
    corrs = np.corrcoef(agg_cstr_dom.reshape(ns, nc*nd).T)

    plt.figure(figsize = [10,10])
    plt.imshow(corrs, aspect='auto', vmin=-1, vmax=1, interpolation='none')

    cstr = ['Det', 'Inc', 'Alt', 'Inf', 'Opn']
    dom  = ['Gen', 'Soc', 'Slf', 'Pol', 'Nws']

    #plt.xticks(np.arange(25), ['{}-{}'.format(c, d) for c in cstr for d in dom], rotation=45, ha='left')
    plt.xticks([])
    plt.yticks(np.arange(25), ['{}-{}'.format(c, d) for c in cstr for d in dom])

    # Insert horizontal lines at every 5th row, vertical at every 5th column
    for i in range(1, 5):
        plt.axhline(y=i*5-0.5, color='k', linewidth=2)
        plt.axvline(x=i*5-0.5, color='k', linewidth=2)

    # Label each cell with numeric value
    nrows, ncols = corrs.shape
    for i in range(nrows):
        for j in range(ncols):
            plt.text(j, i, str(round(corrs[i,j], 1)), ha='center', va='center')

    plt.title('Domain Score Correlation Matrix')
    plt.xlabel('Domain')
    plt.ylabel('Domain')
    #plt.colorbar()
    plt.tight_layout()

def plot_reverse_coding_score_corrs(agg_cstr_dom_rev):
    ns, nc, nd, nr = agg_cstr_dom_rev.shape

    # Flip reverse coded questions
    flipped = agg_cstr_dom_rev.copy()
    #flipped[:,:,:,1] = -agg_cstr_dom_rev[:,:, :, 1]

    corrs = np.corrcoef(flipped.reshape(ns, nc*nd*nr).T)

    cstr = ['Detail', 'Incrementalism', 'Alternatives', 'Info. Seeking', 'Openness']
    dom  = ['Gen', 'Soc', 'Slf', 'Pol', 'Nws']
    rev  = ['Non', 'Rev']
    labels = ['{}-{}-{}'.format(c, d, r) for c in cstr for d in dom for r in rev]

    #plt.figure(figsize = [18,6])
    #plt.suptitle('Reversed Score Correlation Matrix')

    for c in range(5):
        #plt.subplot(1,5,c+1)
        plt.figure(figsize = [4.4,4.4])
        lb, ub = c*10, c*10+10
        plt.imshow(corrs[lb:ub, lb:ub], aspect='auto', vmin=-1, vmax=1, interpolation='none')

        plt.title('{}'.format(cstr[c]))
        plt.xlabel('Domain')
        plt.ylabel('Domain')

        plt.xticks([])
        plt.yticks(np.arange(10), ['{}-{}'.format(d, r) for d in dom for r in rev])

        # Insert horizontal lines at every 2nd row, vertical at every 2nd column
        for i in range(1, 5):
            plt.axhline(y=i*2-0.5, color='k', linewidth=2)
            plt.axvline(x=i*2-0.5, color='k', linewidth=2)

        # Label each cell with numeric value
        nrows, ncols = corrs[lb:ub, lb:ub].shape
        for i in range(nrows):
            for j in range(ncols):
                plt.text(j, i, str(round(corrs[lb:ub, lb:ub][i,j], 1)), ha='center', va='center')

        #plt.colorbar()
        plt.tight_layout()

def plot_domain_score_corrs_separated(agg_cstr_dom):
    ns, nc, nd = agg_cstr_dom.shape

    corrs = np.corrcoef(agg_cstr_dom.reshape(ns, nc*nd).T)

    cstr = ['Detail', 'Incrementalism', 'Alternatives', 'Info. Seeking', 'Openness']
    dom  = ['Gen', 'Soc', 'Slf', 'Pol', 'Nws']
    rev  = ['Non', 'Rev']
    labels = ['{}'.format(d) for d in dom]

    for c in range(5):
        #plt.subplot(1,5,c+1)
        plt.figure(figsize = [4.4,4.4])
        lb, ub = c*5, c*5+5
        plt.imshow(corrs[lb:ub, lb:ub], aspect='auto', vmin=-1, vmax=1, interpolation='none')

        plt.title('{}'.format(cstr[c]))
        plt.xlabel('Domain')
        plt.ylabel('Domain')

        plt.xticks([])
        plt.yticks(np.arange(5), ['{}'.format(d) for d in dom])

        # Insert horizontal lines at every 2nd row, vertical at every 2nd column
        for i in range(1, 5):
            plt.axhline(y=i-0.5, color='k', linewidth=2)
            plt.axvline(x=i-0.5, color='k', linewidth=2)

        # Label each cell with numeric value
        nrows, ncols = corrs[lb:ub, lb:ub].shape
        for i in range(nrows):
            for j in range(ncols):
                plt.text(j, i, str(round(corrs[lb:ub, lb:ub][i,j], 1)), ha='center', va='center')

        #plt.colorbar()
        plt.tight_layout()



# Plot question correlations by question standard deviations
def plot_corrs_by_stds(ansbyq, forward_coded, reverse_coded, keep=None):
    qqcorrs = np.corrcoef(ansbyq)
    avgqcorrs = np.mean(qqcorrs, axis=0)
    qstds = ansbyq.std(axis=1)

    if keep is not None:
        forward_coded = np.intersect1d(forward_coded, keep)
        reverse_coded = np.intersect1d(reverse_coded, keep)

    plt.figure(figsize=[4.67, 4.27])
    plt.plot(qstds[forward_coded], avgqcorrs[forward_coded], 'o')
    plt.plot(qstds[reverse_coded], avgqcorrs[reverse_coded], 'o')
    plt.title('Corrs by SDs')
    plt.xlabel('Standard Deviation')
    plt.ylabel('Average Correlation')
    plt.legend(['Forward Coded', 'Reverse Coded'])
    plt.tight_layout()

# Plot question stds by question means
def plot_stds_by_means(ansbyq, forward_coded, reverse_coded, keep=None):
    qstds = ansbyq.std(axis=1)
    qavgs = np.mean(ansbyq, axis=1)

    if keep is not None:
        forward_coded = np.intersect1d(forward_coded, keep)
        reverse_coded = np.intersect1d(reverse_coded, keep)

    plt.figure(figsize=[4.67, 4.27])
    plt.plot(qstds[forward_coded], qavgs[forward_coded], 'o')
    plt.plot(qstds[reverse_coded], qavgs[reverse_coded], 'o')
    plt.title('Means by SDs')
    plt.xlabel('Standard Deviation')
    plt.ylabel('Average Response')
    plt.legend(['Forward Coded', 'Reverse Coded'])
    plt.tight_layout()

# Plot question correlations by question means
def plot_corrs_by_means(ansbyq, forward_coded, reverse_coded, keep=None):
    qqcorrs = np.corrcoef(ansbyq)
    avgqcorrs = np.mean(qqcorrs, axis=0)
    qavgs = np.mean(ansbyq, axis=1)

    if keep is not None:
        forward_coded = np.intersect1d(forward_coded, keep)
        reverse_coded = np.intersect1d(reverse_coded, keep)

    plt.figure(figsize=[4.67, 4.27])
    plt.plot(qavgs[forward_coded], avgqcorrs[forward_coded], 'o')
    plt.plot(qavgs[reverse_coded], avgqcorrs[reverse_coded], 'o')
    plt.title('Correlations by Means')
    plt.xlabel('Average Response')
    plt.ylabel('Average Correlation')
    plt.legend(['Forward Coded', 'Reverse Coded'])
    plt.tight_layout()


def plot_qstds(qstats, forward_coded, reverse_coded, keep=None):

    sorted_inds = np.argsort(qstats['std'])
    if keep is not None:
        sorted_inds = np.array([i for i in sorted_inds if i in keep])

    rc_msk = [i in reverse_coded for i in sorted_inds]
    fc_msk = [i in forward_coded for i in sorted_inds]
    rc_ind = np.where(rc_msk)[0]
    fc_ind = np.where(fc_msk)[0]

    std_half = 0.82
    std_unif = 1.414
    std_ends = 2.0

    plt.figure(figsize=[4.67, 4.27])
    plt.plot(fc_ind, qstats['std'].values[sorted_inds[fc_msk]], 'o')
    plt.plot(rc_ind, qstats['std'].values[sorted_inds[rc_msk]], 'o')
    plt.axhline(std_half, color='r', linestyle='--', label='Half')
    plt.axhline(std_unif, color='g', linestyle='--', label='Uniform')
    plt.axhline(std_ends, color='b', linestyle='--', label='Ends')
    plt.title('Question SDs')
    plt.xlabel('Question')
    plt.ylabel('Standard Deviation')
    plt.grid()
    plt.legend()
    plt.tight_layout()


def plot_running_std(running_std):
    std_half = 0.82
    std_unif = 1.414
    std_ends = 2.0

    grand_std = np.mean(running_std, axis=0)
    plt.figure(figsize=[4.67, 4.27])
    plt.axhline(std_half, color='r', linestyle='--', label='Half')
    plt.axhline(std_unif, color='g', linestyle='--', label='Uniform')
    plt.axhline(std_ends, color='b', linestyle='--', label='Ends')
    plt.legend()
    plt.plot(grand_std, '--k', linewidth=2)
    plt.plot(running_std.T, alpha=0.2)
    plt.title('Running Standard Deviations')
    plt.xlabel('Question Group')
    plt.ylabel('Standard Deviation')
    plt.xlim(0, 91)
    plt.tight_layout()


def plot_running_average_trends(running_avg, sids):
    import statsmodels.api as sm
    subject_models = []
    params = np.zeros((len(sids), 2))
    pvals  = np.zeros((len(sids), 1))
    for i in range(len(sids)):
        x = np.arange(0, 91)
        y = running_avg[i,:]
        x = sm.add_constant(x)
        model = sm.OLS(y, x).fit()
        params[i,0] = model.params[0]
        params[i,1] = model.params[1]
        pvals[i] = model.pvalues[1]


    # Plot distribution of slopes
    plt.figure(figsize=[4.67, 4.27])
    plt.hist(params[:,1]*91, bins=20)
    plt.title('Distribution of SD by Time Slopes')
    plt.xlabel('Slope')
    plt.ylabel('Count')
    plt.xticks(rotation=45)
    plt.axvline(0, color='r', linestyle='--', label='Zero Slope')
    plt.tight_layout()

    # Plot distribution of p-values
    plt.figure(figsize=[4.67, 4.27])
    plt.hist(pvals, bins=np.arange(0, 1.01, 0.01), label='p-values')
    plt.title('Distribution of p-values')
    plt.xlabel('p-value')
    plt.ylabel('Count')
    plt.xticks(rotation=45)
    plt.axvline(0.05, color='r', linestyle='--', label='p=0.05')
    plt.legend()
    plt.tight_layout()

def plot_running_std_trends(running_std, sids):
    # Get the trend over time in standard deviations for each subject using a linear model
    import statsmodels.api as sm
    subject_models = []
    params = np.zeros((len(sids), 2))
    pvals  = np.zeros((len(sids), 1))
    for i in range(len(sids)):
        x = np.arange(0, 91)
        y = running_std[i,:]
        x = sm.add_constant(x)
        model = sm.OLS(y, x).fit()
        params[i,0] = model.params[0]
        params[i,1] = model.params[1]
        pvals[i] = model.pvalues[1]

    # Plot distribution of SD slopes
    plt.figure(figsize=[4.67, 4.27])
    plt.hist(params[:,1]*91, bins=20)
    plt.title('Distribution of SD by Time Slopes')
    plt.xlabel('Slope')
    plt.ylabel('Count')
    plt.xticks(rotation=45)
    plt.axvline(0, color='r', linestyle='--', label='Zero Slope')
    plt.tight_layout()

    # Plot distribution of p-values
    plt.figure(figsize=[4.67, 4.27])
    plt.hist(pvals, bins=np.arange(0, 1.01, 0.01), label='p-values')
    plt.title('Distribution of p-values')
    plt.xlabel('p-value')
    plt.ylabel('Count')
    plt.xticks(rotation=45)
    plt.axvline(0.05, color='r', linestyle='--', label='p=0.05')
    plt.legend()
    plt.tight_layout()



# Plot running averages
def plot_running_average_responses(running_avg):
    grand_avg = np.mean(running_avg, axis=0)
    plt.figure(figsize=[4.67, 4.27])
    plt.plot(grand_avg, '--', color='k', label='Grand Average', linewidth=2)
    plt.plot(running_avg.T, alpha=0.2)
    plt.title('Running Averages')
    plt.xlabel('Question Group')
    plt.ylabel('Average')
    plt.xlim(0, 91)
    plt.tight_layout()

def plot_participant_grand_stds(ansbyq):

    # Get average question standard deviation by participant
    participant_std = np.std(ansbyq, axis=0)
    sorted_participant_std = np.sort(participant_std)

    std_half = 0.82
    std_unif = 1.414
    std_ends = 2.0

    plt.figure(figsize=[4.67, 4.27])
    plt.plot(sorted_participant_std, '-o')
    plt.axhline(std_half, color='r', linestyle='--', label='Half')
    plt.axhline(std_unif, color='g', linestyle='--', label='Uniform')
    plt.axhline(std_ends, color='b', linestyle='--', label='Ends')
    plt.title('Participant Grand SDs')
    plt.xlabel('Participant')
    plt.ylabel('Standard Deviation')
    plt.grid()
    plt.legend()
    plt.tight_layout()


def plot_construct_PCAs(constructs, qcoords, zbyq):
    ve = []
    top_pcs = []
    for i in range(5):
        inds = np.where(qcoords['cnum'] == i)[0]
        pca = PCA()
        pca.fit(zbyq[inds,:].T)
        ve.append(pca.explained_variance_ratio_[0:10])
        top_pcs.append(pca.components_[0,:])

    ve = np.array(ve)
    top_pcs = np.array(top_pcs)

    # plt.figure(figsize=[4.67, 4.27])
    # plt.plot(ve[:,0], '-o')
    # plt.title('VE by Top PC for Each Construct')
    # plt.xlabel('Construct')
    # plt.ylabel('Fractional Variance Explained')
    # plt.xticks(np.arange(0,5), constructs)
    # plt.grid('on')
    # plt.tight_layout()

    # Plot within-construct variance explained
    plt.figure(figsize=[4.67, 4.27])
    plt.plot(ve.T, '-o')
    plt.title('Within-Construct VEs')
    plt.legend(constructs)
    plt.xlabel('PC')
    plt.ylabel('Fractional Variance Explained')
    plt.tight_layout()

    # Plot within-construct PC loadings
    plt.figure(figsize=[4.67, 4.27])
    plt.plot(top_pcs.T, '-o')
    plt.title('Top PC Loadings by Construct')
    plt.plot([0,20], [0,0], 'k--')
    plt.legend(constructs)
    plt.xlabel('Question')
    plt.ylabel('Loading')
    plt.tight_layout()

def plot_construct_score_stds(agg_cstr, constructs):
    plt.figure(figsize=[4.67, 4.27])
    plt.plot(np.std(agg_cstr,axis=0), '-o')
    plt.title('Construct Score Standard Deviations')
    plt.xlabel('Construct')
    plt.ylabel('Standard Deviation')
    #plt.grid('on')
    plt.xticks(np.arange(0,5), constructs)
    plt.tight_layout()

# ---- Plot PCA of domain scores ---- #
def plot_domain_score_pca(agg_cstr_dom, qcoords):
    sbyd = agg_cstr_dom.reshape(76, 25)
    pca = PCA()
    pca.fit(sbyd)
    for i in range(3):
        pc = pca.components_[i].reshape(5,5)
        plot_construct_domain_heatmap(pc, qcoords, title='PC{} Loadings'.format(i+1))

    # Plot domain-score variance explained
    plt.figure(figsize=[4.67, 4.27])
    plt.plot(pca.explained_variance_ratio_, '-o')
    plt.title('PCA of Domain Aggregated Scores')
    plt.xlabel('PC')
    plt.ylabel('Fractional Variance Explained')
    plt.tight_layout()


def plot_component_heatmaps(pca, qcoords, ncomps):
    pca_components = pca.components_[0:ncomps,:].T
    for i in range(ncomps):
        pc_heatmap = create_construct_domain_item_heatmap(pca_components[:,i], qcoords)
        plot_construct_domain_item_heatmap(pc_heatmap, title='PC{} Loadings'.format(i))



def plot_pca_of_subsetted_questions(ansbyq, qcoords, drop_idx, ncomps=3):
    remaining = np.setdiff1d(np.arange(100), drop_idx)
    reduced_abq = ansbyq[remaining,:]

    # Get new PCA of the data
    pca = PCA()
    pca.fit(reduced_abq.T)

    plot_pca_variance_explained(pca, title='Item Level PCA Variance Explained')
    plot_pca_components(pca, xlabel='Item', title='Item Level PCA Components')

    pca_components = pca.components_[0:ncomps,:].T
    for i in range(ncomps):
        pc = nanfill_to_original(pca.components_[i], remaining)
        pc_heatmap = create_construct_domain_item_heatmap(pc, qcoords)
        plot_construct_domain_item_heatmap(pc_heatmap, title='PC{} Loadings'.format(i))


def plot_qq_corrs(ansbyq, keep=None):
    if keep is not None:
        ansbyq = ansbyq[keep,:]
    
    qqcorrs = np.corrcoef(ansbyq)
    plt.figure()
    plt.imshow(qqcorrs, aspect='auto')
    plt.colorbar()
    plt.xlabel('Item')
    plt.ylabel('Item')
    plt.title('Item-Item Correlation Matrix')