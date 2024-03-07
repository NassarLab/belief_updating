import os
import re
import json
import pandas as pd
from pathlib import Path

# Custom error for non-directory inputs
class NotADirectoryError(Exception):
    pass

# Read demographic information
def read_demographic_info(filename):
    return pd.read_csv(filename)

# Match and append demographic information to participant data
def append_demographic_info(q_processed, demographic_df, prolific_pid):
    participant_info = demographic_df[demographic_df['Participant id'] == prolific_pid]
    if not participant_info.empty:
        participant_info_row = participant_info.iloc[0]  # Assuming each PID has only one row of demographic info
        for q in q_processed:  # Iterate through each question processed
            for col in participant_info.columns:
                if col != 'Participant id':  # Avoid duplicating the participant id
                    q['Demographic_' + col] = participant_info_row[col]
    return q_processed


# Instruction strings for data cleaning
start_str_1 = "<p style='font-size:36px; line-height:1.5; text-align:center;'>Welcome to this part of the experiment! Please note:</p><p style='font-size:36px; line-height:1.5; text-align:left;'><ol style='font-size:36px; line-height:1.5; text-align:left;'><li>It will take about ~40 minutes to complete.</li><li>Please do <b>NOT</b> close the window, refresh the page, or use the back button, as that may terminate the experiment.</li><li>Before you start, please make sure your internet connection is stable, which is very important for successfully completing the experiment.</li></ol></p>"
start_str_2 = "<p style='font-size:36px; line-height:1.5;'>In this part, you will first complete 2 short tasks.<br>Then, you will be asked to fill out 3 questionnaires.<br>Please press the button below to start.</p>"

# The root directory containing the sub-directories
rootdir = "./jresults"

def read_metadata(metadata_path):
    with open(metadata_path, 'r') as file:
        metadata = json.load(file)
    prolific_pid_map = {str(result['id']): result['urlQueryParameters']['PROLIFIC_PID'] for data in metadata['data'] for result in data['studyResults']}
    print(prolific_pid_map)  # Debug print to verify the map
    return prolific_pid_map

def scan_dir(rootdir, metadata_path, demographic_file):
    prolific_pid_map = read_metadata(metadata_path)  # Read metadata to get PROLIFIC_PID map
    demographic_df = read_demographic_info(demographic_file)  # Read demographic information
    num_participants = 0

    for path in Path(rootdir).iterdir():
        if path.is_dir():
            num_participants += 1
            pData = ""
            pID = str(path).split('/')[-1].replace('study_result_', '')  # Extract participant ID from path

            # Get PROLIFIC_PID for current participant
            prolific_pid = prolific_pid_map.get(pID, "Unknown")
            print(f"Processing participant {pID} with Prolific ID: {prolific_pid}")

            trails = os.listdir(path)
            for t in trails:
                data_path = os.path.join(str(path), t, "data.txt")
                if os.path.exists(data_path):
                    with open(data_path, 'r') as infile:
                        data = infile.read()

                        # Process Data
                        cleaned_data = data.replace(start_str_1, '').replace('}', '}\n').replace(start_str_2, '')
                        cleaned_data = re.sub(r'[{}\[\]\"]', '', cleaned_data)
                        pData += cleaned_data + "\n"

            pData = "," + pData[5:]
            lines = pData.split("\n")
            search_string = "qns_items:"
            questionnaire_head = next((i for i, line in enumerate(lines, 1) if line.startswith('item:')), None)
            end_string = next((i for i, line in enumerate(lines, 1) if line.startswith(search_string)), len(lines))

            lines_to_delete = set(range(0, questionnaire_head)) | set(range(end_string, len(lines)))
            lines = [line for index, line in enumerate(lines, start=1) if index not in lines_to_delete]

            q_processed = process_lines(lines, prolific_pid)  # Process lines with prolific_pid
            q_processed = append_demographic_info(q_processed, demographic_df, prolific_pid)  # Append demographic info

            # Write generalization data
            directory_name = "response"
            os.makedirs(directory_name, exist_ok=True)

            df = pd.DataFrame(q_processed)
            file_path = os.path.join(directory_name, f'p{pID}.csv')
            df.to_csv(file_path, index=False)

    print(f"Processed {num_participants} participants.")


def find_first_line_number(strings_list, search_string):
    for line_number, line in enumerate(strings_list, 1):
        if line.startswith(search_string):
          # print(f"'{line_number}' is the first line")
          return line_number
        else:
          continue

def find_line_number(strings_list, search_string):
    for line_number, line in enumerate(strings_list, 1):
        if line.startswith(search_string):
          # print(f"'{search_string}' found on line {line_number}")
          return line_number
        else:
          continue

# Read metadata.json to fetch PROLIFIC_PID for each participant
def read_metadata(metadata_path):
    with open(metadata_path, 'r') as file:
        metadata = json.load(file)
    prolific_pid_map = {str(result['id']): result['urlQueryParameters']['PROLIFIC_PID'] for data in metadata['data'] for result in data['studyResults']}
    return prolific_pid_map

def process_lines(lines, prolific_pid):
    q_processed = []
    for x in lines:
        row_data = {"PROLIFIC_PID": prolific_pid}  # Add the PROLIFIC_PID to each row
        if '</b>' in x:
            end_q_index = x.index('</b>')
            question = x[5:end_q_index]  # Assume there's only one question per line
            row_data["Question"] = question
            items = x[end_q_index+5:].split(',')
            for item in items:
                data = item.split(":")
                if len(data) == 2:
                    row_data[data[0].strip()] = data[1].strip()
            q_processed.append(row_data)
        else:
            items = x.split(',')
            for item in items:
                data = item.split(":")
                if len(data) == 2:
                    row_data[data[0].strip()] = data[1].strip()
            q_processed.append(row_data)
    return q_processed




# Usage example
metadata_path = "./metadata.json"  # Path to your metadata.json file
demographic_file = 'demographic_information.csv'  # Path to your demographic_information.csv file
scan_dir(rootdir, metadata_path, demographic_file)

  # q_processed = []
      # for x in lines[1:]:
      #   if '</b>' in x:
      #     end_q_index = x.index('</b>')
      #     q_processed.append(x[5:end_q_index+4].split('</b>') + x[end_q_index+5:].split(','))
      #   else:
      #     q_processed.append(x.split(','))
      #     df = pd.DataFrame([x.split(',') for x in lines])

      # df = pd.DataFrame(q_processed)
      # df.to_csv(f'Participant_{pID}.csv', index=False)

