import os
import pandas as pd
import numpy as np

def read_response_data(dirstr):
    """
    Cycles over csv files in a directory, reading subject data into a stacked DF.
    """
    # Get data files
    files = sorted(os.listdir(dirstr))
    files = [file for file in files if file.startswith('p') and file.endswith('.csv')]

    # Initialize drop counter
    drop_cnt = 0

    # Manually marked as bad list
    #bad_list = list(pd.read_csv(dir + '/remove_list.csv', header=1))
    bad_list = [13619, 13624, 13570] #[13407]
    drop_cnt += len(bad_list)
    snumlen = 5

    expected_nq = 104

    # Read all the data
    sids, data = [], []
    print('Reading in data from ' + dirstr)
    for file in files:
        
        # Get subject number from filename
        #sid = int(file[1:snumlen+1])
        sid = int(file[5:snumlen+5])

        # Check if bad
        if sid in bad_list: continue

        # Data from file, append sid
        df = pd.read_csv(dirstr + file)
        df['sid'] = sid

        # Check that all subjects have same # of questions
        nq = df.shape[0]
        if nq != expected_nq:
            print('Subject '+ str(sid) + ' has ' + str() + ' questions.')

        # Check attention failures
        failures = attention_failures(df)

        # Notify
        if failures > 0:
            print('Subject ' + str(sid) + ' failed ' + str(failures) + ' attention checks, dropping.')
            drop_cnt +=1
            continue

        # Save to list
        sids.append(sid)
        data.append(df)

    # Notify total # of dropped subjects
    print('Dropped ' + str(drop_cnt) + ' of ' + str(len(files)) + ' subjects.')
    print('Fraction remaining: ' + str((len(files) - drop_cnt)/len(files)))

    # Merge data into single frame
    data = pd.concat(data).reset_index(drop = True)

    # Make sure we filtered properly
    assert not any([sid in bad_list for sid in sids])

    # Convert question number to an actual number
    data.quest_num = data.quest_num.apply(lambda x: int(x.split('/')[0]))

    return sids, data


def attention_failures(df):
    """
    Quality controls for the old data, needs updating for the new.
    """
    # Attention check questions and answers
    checks  = [15,35,55,75]
    answers = [ 1, 1, 3, 4]

    # Verify that these are as expected
    attn_qs = np.where(df['is_attention'])[0]
    assert np.all(attn_qs == checks)

    # Check if subject passed them
    failures = 0
    for i, check in enumerate(checks):
        if df.iloc[check].answer_num != answers[i]:
            failures += 1

    return failures