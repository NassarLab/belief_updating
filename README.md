# Belief Updating Questionnaire Analysis

This repository contains files for analyzing belief updating questionnaire responses.

### To run the code: ###
- Put subject .csv files in a local folder ```./data/```
- If this is the take-2 pilot data, put it under ```./data/take-2-pilot/```
- If not, change the folder name in configs.py.
- Open ipython in the terminal (or use ```pip install ipython```, then do so)
- Use ```%run run.py``` to run all the currently implemented analyses

### File contents: ###
- **analysis.py**: contains the main data analysis tools
- **configs.py**: imports and configuration for using this package
- **depr.py**: deprecated code, snippets that
- **plots.py**: functions for plotting results of analyses
- **readin.py**: functions for reading subject data and performing QCs
- **run.py**: script for reading subject data, running analyses, and plotting things
- **utils.py**: misc. functions