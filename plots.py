import numpy as np
import matplotlib.pyplot as plt

def plot_response_dists(dfs, inds, title, short_qs):
    # Plot response DFS
    fig, ax = plt.subplots(figsize = [8,10])
    cax = ax.matshow(dfs[inds,:], aspect='auto', vmin=0, vmax=1)

    plt.title(title)

    # Set questions as y-axis tick labels
    ax.yaxis.set_ticks_position('right')
    ax.set_yticks(inds-inds[0])
    ax.set_yticklabels([short_qs[i] for i in inds])

    # Set agreement as x-axis tick labels
    ax.set_xticks([0,1,2,3])
    ax.set_xticklabels(['Strongly Disagree','Slightly Disagree', 'Slightly Agree', 'Strongly Agree'], rotation = 30, ha='left')

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

def plot_pca_variance_explained(pca):
    # Plot relative variance explained
    plt.figure(figsize = [4,4])
    plt.plot(pca.explained_variance_ratio_, 'o')
    plt.title('PCA Scree Plot')
    plt.xlabel('Component')
    plt.ylabel('Variance Explained')
    plt.grid('on')
    plt.tight_layout()

def plot_pca_cumulative_variance(pca):
    # Plot cumulative variance explained
    plt.figure(figsize = [4,4])
    plt.plot(np.cumsum(pca.explained_variance_ratio_),'o')
    plt.title('PCA Cumulative Fractions')
    plt.xlabel('Component')
    plt.ylabel('Cumulative Relative Variance')
    plt.grid('on')
    plt.tight_layout()

def plot_pca_components(pca, comps=[0,1]):
    # Plot cumulative variance explained
    plt.figure(figsize = [4,4])
    plt.plot(pca.components_[comps,:].T,'-o')
    plt.title('PCA Components')
    plt.xlabel('Question')
    plt.ylabel('Loading')
    plt.grid('on')
    plt.tight_layout()

def plot_corrs(data):
    corrs = np.corrcoef(data)
    plt.matshow(corrs, aspect='auto')
    plt.title('Answer Correlation Matrix')
    plt.xlabel('Question')
    plt.ylabel('Question')
    plt.colorbar()
    plt.tight_layout()

def plot_construct_domain_heatmap(heatmap):
    """
    Plots a heatmap on the excel grid. Heatmap should be on this grid already.
    """
    # Plot the data
    plt.matshow(heatmap, aspect='auto')
    plt.title('Construct-Domain Heatmap')

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
