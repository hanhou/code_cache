# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 19:53:41 2021

Han based on Marton's code
https://github.com/ilyakolb/2P_cellAttached_pipeline/blob/master/notebook_google/notebook_main.py

Use Offce365-REST-Python-Client
https://github.com/vgrem/Office365-REST-Python-Client
https://stackoverflow.com/questions/48424045/how-to-read-sharepoint-online-office365-excel-files-in-python-with-work-or-sch

*** However, this does not work due to Multi-factor Authentication (MFA) constraint of Allen's sharepoint settings
*** Workaround: Use the "sync" functionality of Teams to sync the metadata to the local PC (much easier!)
"""

import json
import pandas as pd
import time
import os.path


# ======== Settings =======
# Settings

local_sharepoint_sync_folder = R'C:/Users/admin/Allen Institute/Neural Dynamics - metadata'
remote_meta_file_animal = R'Surgery, water restriction and training.xlsx'
remote_meta_file_lab = R'Lab metadata.xlsx'
    
# Target folder
with open('D:/Han_Sync/Svoboda/Scripts/map-ephys/dj_local_conf.json') as f:
        dj_conf = json.load(f)
local_meta_dir_animal = dj_conf['custom']['behavior_bpod']['meta_dir']
local_meta_dir_lab = dj_conf['custom']['behavior_bpod']['meta_lab_dir']
    
   
def remote_meta_to_csvs(remote_notebook_name, local_meta_dir, transposed = False):
    # Save remote notebook to local csv files, one for each sheet
    # Update anyway (gspread with v4 api is not easy to get the last modified time)

    excel_file = pd.ExcelFile(os.path.join(local_sharepoint_sync_folder, remote_notebook_name))
    sheet_names = excel_file.sheet_names

    for sheet in sheet_names:
        df_wr = pd.read_excel(excel_file, sheet)
                
        if type(df_wr) == pd.DataFrame:
            df_wr.to_csv(os.path.join(local_meta_dir, f'{sheet}.csv')) 
            print(f'  csv saved: {sheet} in {remote_notebook_name}')
        else:
            print(f'  csv failed to save: {sheet} in {remote_notebook_name}')
        
        
def fetch_foraging_meta_from_sharepoint():

    print('Updating surgery and WR metadata from Sharepoint sync folder...')
    remote_meta_to_csvs(remote_meta_file_animal, local_meta_dir_animal)
    print('done!')
        
    print('Updating Lab metadata from Sharepoint sync folder...')
    remote_meta_to_csvs(remote_meta_file_lab, local_meta_dir_lab)
    print('done!')
 
if __name__ == '__main__':
    fetch_foraging_meta_from_sharepoint()