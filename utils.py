import numpy as np
import pandas as pd

def flip_reverse_coded(ansbyq, reverse_coded):
    """
    Flips the answers of reverse-coded questions.
    """
    recoded = np.copy(ansbyq)
    # Flip reverse coded questions
    for qidx in reverse_coded:
        recoded[qidx,:] = 3 + -(ansbyq[qidx,:] - 3)

    return recoded


def create_construct_domain_heatmap(data):
    """
    Gets coordinates into our excel grid and creates a heatmap on this grid.
    """

    # Heatmap of the data, reorganized into 5 constructs with 5 domains and 4 questions each
    heatmap = np.zeros([5*4, 5])

    qcoords = np.zeros([100, 5])
    for qnum in range(100):
        # Get construct, domain, and row number for grid
        cnum = np.floor(qnum /  20).astype(int)
        dnum = np.floor(qnum % 20 / 4).astype(int)
        rnum = (qnum % 20) % 4

        # Get coordinates into full grid
        i, j = cnum*4 + rnum, dnum

        # Used for debugging
        #print(qnum, cnum, dnum, rnum, i, j)

        # Regrid the data, save the q-coordinates
        heatmap[i,j] = data[qnum]
        qcoords[qnum, :] = (cnum, dnum, rnum, i, j)

    # Save labels for qcoords
    qcoords = pd.DataFrame(qcoords, columns = ['cnum', 'dnum', 'rnum', 'i', 'j'], dtype = int)

    return heatmap, qcoords