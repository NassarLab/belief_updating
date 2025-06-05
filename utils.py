import numpy as np
import pandas as pd
import scipy as sp

def flip_reverse_coded(ansbyq, reverse_coded):
    """
    Flips the answers of reverse-coded questions.
    """
    recoded = np.copy(ansbyq)
    # Flip reverse coded questions
    for qidx in reverse_coded:
        recoded[qidx,:] = 3 + -(ansbyq[qidx,:] - 3)

    return recoded

def get_general_question_coordinates():
    """
    Create maps from question numbers to...
    1) construct, domain, and reverse-coding (i.e. hierarchical levels)
    2) Row number and grid indices in construct-domain grid
    """
    # Manually generated list of reverse coded questions
    reverse_coded = np.array([
        1 ,3 , 5 ,6 , 11,12, 14,15, 18,20,
        22,23, 27,28, 30,31, 35,36, 38,40,
        41,44, 46,48, 51,52, 54,55, 57,59,
        62,64, 67,68, 71,72, 75,76, 78,80,
        83,84, 86,87, 90,92 ,94,96, 98,100]) -1

    qcoords = np.zeros([100, 6])
    for qnum in range(100):
        # Get construct, domain, and row number for grid
        cnum = np.floor(qnum /  20).astype(int)
        dnum = np.floor(qnum % 20 / 4).astype(int)
        rnum = (qnum % 20) % 4

        # Get coordinates into full grid
        i, j = cnum*4 + rnum, dnum

        # Used for debugging
        #print(qnum, cnum, dnum, rnum, i, j)

        # Save the q-coordinates
        rev = 1 if qnum in reverse_coded else 0
        qcoords[qnum, :] = (cnum, dnum, rev, rnum, i, j)
        
    # Save labels for qcoords
    qcoords = pd.DataFrame(qcoords, columns = ['cnum', 'dnum', 'rev', 'rnum', 'i', 'j'], dtype = int)

    return qcoords


def create_construct_domain_item_heatmap(data, qcoords):
    """
    Gets coordinates into our excel grid and creates a heatmap on this grid.
    """
    # Heatmap of the data, reorganized into 5 constructs with 5 domains and 4 questions each
    heatmap = np.zeros([5*4, 5])
    for qnum in range(100):
        i,j = qcoords.iloc[qnum][['i', 'j']]
        # Regrid the data, save the q-coordinates
        heatmap[i,j] = data[qnum]

    return heatmap


def ans_to_scores(ansbyq, qcoords, flip=False):
    """
    Converts the answer matrix to scores by domain and construct.

    If reverse coded are already flipped, flip should be false.
    """
    # Z-scores for each question
    zbyq = sp.stats.zscore(ansbyq, axis=1)
    ns   = zbyq.shape[1]

    # Initialize aggregate scores over each domain, construct, reversal
    agg_cstr_dom_rev = np.zeros([ns, 5, 5, 2])
    agg_cstr_dom     = np.zeros([ns, 5, 5])
    agg_cstr         = np.zeros([ns, 5])

    # For each question, get the coordinates and add to aggregates
    for qnum in range(100):

        # Get the coordinates and scores for this question
        cnum, dnum, rev = qcoords.iloc[qnum][['cnum', 'dnum', 'rev']]
        scores = zbyq[qnum,:]

        # Flip the reverse-coded scores
        if rev == 1 and flip:
            scores = -scores
            zbyq[qnum,:] = scores

        # Add to aggregates
        agg_cstr[:, cnum] += scores
        agg_cstr_dom[:, cnum, dnum] += scores
        agg_cstr_dom_rev[:, cnum, dnum, rev] += scores

    return zbyq, agg_cstr, agg_cstr_dom, agg_cstr_dom_rev

def get_cdf(data, type = 'cdf'):
    vals = np.unique(data)
    cnts = np.array([sum(item == data) for item in vals])
    csum = np.cumsum(cnts)
    cdf  = csum/csum[-1]
    pdf  = cnts/sum(cnts)

    df = cdf if type == 'cdf' else pdf

    return vals, df

def nanfill_to_original(vec, remaining):
    new_vec = np.full(100, np.nan)
    new_vec[remaining] = vec
    return new_vec