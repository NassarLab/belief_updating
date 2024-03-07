
def do_round_analysis(round):

    # Get subject response data
    sids, data = read_response_data()

    # Get subject attention check performance
    sids, data = quality_check(sids, data)

    # Get question data
    qs, ansbyq, qstats, qcdfs, qpdfs, inds = get_qdata(data)

    # Compare baseline and reduced questionnaire PCs
    keep = run_pca_comparison(ansbyq)

    # Plot answer means and standard deviations
    plot_means(qstats)
    plot_stds(qstats)

    # Plot PCA results 
    plot_pca_variance_explained( pca)
    plot_pca_cumulative_variance(pca)

    # Get abbreviated questions for plots
    short_qs = [' '.join(Q.split(' ')[0:10]) for Q in qs]

    # Plot top and bottom 40 most informative PDFs
    plot_response_dists(qpdfs, np.arange(0,40)   , 'Response PDFs (Top 40)', short_qs)
    plot_response_dists(qpdfs, np.arange(110,150), 'Response PDFs (Bottom 40)', short_qs)