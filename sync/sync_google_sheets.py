# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 19:53:41 2021

Han based on Marton's code
https://github.com/ilyakolb/2P_cellAttached_pipeline/blob/master/notebook_google/notebook_main.py
"""

import gspread
import json
import pandas as pd
import time
import os.path

is_windows = 'nt' in os.name

# ======== Settings =======
# Remote
if is_windows:
    gspread_client = gspread.service_account(filename='F:/Data_for_ingestion/gspread_key.json')
else:
    gspread_client = gspread.service_account(filename='/mnt/f/Data_for_ingestion/gspread_key.json')
    
remote_meta_file_animal = 'Surgery, water restriction and training'
remote_meta_file_lab = 'Lab metadata'

# Local
if is_windows:
    with open('D:/Han_Sync/Svoboda/Scripts/map-ephys/dj_local_conf.json') as f:
        dj_conf = json.load(f)
    local_meta_dir_animal = dj_conf['custom']['behavior_bpod']['meta_dir']
    local_meta_dir_lab = dj_conf['custom']['behavior_bpod']['meta_lab_dir']
else:
    local_meta_dir_animal = '/mnt/f/Data_for_ingestion/Foraging_behavior/Metadata'
    local_meta_dir_lab = '/mnt/f/Data_for_ingestion/Foraging_behavior/Metadata_lab'
    

def fetch_sheet_titles(notebook_name):
    notebook = gspread_client.open(notebook_name)
    worksheets = notebook.worksheets()
    return [sheet.title for sheet in worksheets]

def fetch_sheet(notebook_name, sheet_title, transposed = False):
    notebook = gspread_client.open(notebook_name)
        
    data = notebook.worksheet(sheet_title).get_all_records()
    try:
        df = pd.DataFrame(data)
        return df
    except:
        print([notebook_name,sheet_title])
        print(data)
        return None
    
def remote_meta_to_csvs(remote_notebook_name, local_meta_dir, transposed = False):
    # Save remote notebook to local csv files, one for each sheet
    # Update anyway (gspread with v4 api is not easy to get the last modified time)
    sheets = fetch_sheet_titles(remote_notebook_name)
    for sheet in sheets:
        while True:
            try:    
                df_wr = fetch_sheet(remote_notebook_name, sheet, transposed)
                break
            except gspread.exceptions.APIError as err:
                print(err)
                print(f'quota exceeded at table {sheet}, waiting 150s')
                time.sleep(150)
                
        if type(df_wr) == pd.DataFrame:
            df_wr.to_csv(os.path.join(local_meta_dir,f'{sheet}.csv')) 
            print(f'  csv saved: {sheet} in {remote_notebook_name}')
        else:
            print(f'  csv failed to save: {sheet} in {remote_notebook_name}')
        
        
def fetch_foraging_meta_from_goolge():

    print('Updating surgery and WR metadata from google drive...')
    remote_meta_to_csvs(remote_meta_file_animal, local_meta_dir_animal)
    print('done!')
        
    print('Updating Lab metadata from google drive...')
    remote_meta_to_csvs(remote_meta_file_lab, local_meta_dir_lab)
    print('done!')
 
if __name__ == '__main__':
    fetch_foraging_meta_from_goolge()