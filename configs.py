# Imports
import os
import sklearn as sk
import numpy   as np
import pandas  as pd
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA, SparsePCA

# Package functions
from readin   import *
from analysis import *
from plots    import *
from utils    import *

# Interactive plotting
plt.ion()

# Subject data folder
subj_data_dir = './data/take-2-pilot/'