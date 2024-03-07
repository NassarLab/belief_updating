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