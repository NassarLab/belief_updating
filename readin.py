import os
import pandas as pd

def import_data(dirstr, round = 3):
    """
    Wrapper for importing data then doing QCs. 
    QCs aren't implemented for new data yet.
    """
    # Read subject data
    sids, data = read_response_data(dirstr)
    
    # Get subject attention check performance
    if round in [1,2]:
        sids, data = quality_check(sids, data, round)

    return sids, data

def read_response_data(dirstr):
    """
    Cycles over csv files in a directory, reading subject data into a stacked DF.
    """
    # Get data files
    files = sorted(os.listdir(dirstr))
    files = [file for file in files if file.startswith('p') and file.endswith('.csv')]

    # Manually marked as bad list
    #bad_list = list(pd.read_csv(dir + '/remove_list.csv', header=1))
    bad_list = []
    snumlen = 5

    # Read all the data
    sids, data = [], []
    print('Reading in data from ' + dirstr)
    for file in files:

        # Check if bad
        if file[1:snumlen+1] in bad_list: continue

        # Subject id from filename
        sids.append(file[1:snumlen+1])

        # Data from file, append sid
        df = pd.read_csv(dirstr + file)
        df['sid'] = sids[-1]

        # The validation set (weirdly) has an extra NAN row.
        df = df.dropna(subset=['Question'])

        # Insert basic check that all subjects have same # of questions
        print('Subject '+ sids[-1] + ' has ' + str(df.shape[0]) + ' questions.')

        # Save to list
        data.append(df)

    # Merge data into single frame
    data = pd.concat(data).reset_index(drop = True)

    # Make sure we filtered properly
    assert not any([sid in bad_list for sid in sids])

    # Convert question number to an actual number
    data.quest_num = data.quest_num.apply(lambda x: int(x.split('/')[0]))

    return sids, data


def quality_check(sids, data, round):
    """
    Quality controls for the old data, needs updating for the new.
    """
    # Three check questions
    cqn = [22,71,115] if round == 1 else [13,42,71]
    check_1 = list(data.loc[data.quest_num == cqn[0],:].answer_num == 3)
    check_2 = list(data.loc[data.quest_num == cqn[1],:].answer_num == 1)
    check_3 = list(data.loc[data.quest_num == cqn[2],:].answer_num == 1)

    # Get run lengths for each subject
    runlens, check_4 = [], []
    for i, s in enumerate(sids):

        # Detect changes in answer
        df = pd.DataFrame()
        df['shifted'] = data[data.sid == s]['answer_num'].shift(1) != data[data.sid == s]['answer_num']

        # Cumulative sum of bools tells us which chunk (run) each answer falls in
        df['chunk'] = df['shifted'].cumsum()

        # Group them by run and count how many are in each run
        runlens.append( df.groupby('chunk').size().tolist() )

        # Check if any runs are longer than 10
        check_4.append( any([l < 10 for l in runlens[i]]) )

    # TODO: Add variance check back in?

    # List of pass/fail for each subject
    pass_check = []
    for i in range(len(check_1)):
        pass_check.append(check_1[i] and check_2[i] and check_3[i] and check_4[i])

    # Failed subject list
    failed = [sids[i] for i, val in enumerate(pass_check) if not val]

    # Display who failed
    if len(failed) == 0:
        print('No subjects failed attention checks.')
    else:
        print('Subjects failing checks:')
        print(failed)

    # Remove subjects failing checks from data
    for sid in failed:
        inds = data.index[data.sid == sid].tolist()
        data = data.drop(inds)
        data = data.reset_index(drop = True)

    # Remove subjects failing checks from sids list
    sids = [sid for sid in sids if sid not in failed]

    return sids, data