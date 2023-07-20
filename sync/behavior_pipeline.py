'''
Automatic pipeline for behavioral ingestion
1. Copy bpod sessions from all the rigs
2. Sync spreadsheet
3. Update lab meta (mice, person, ...)
4. Ingest behavioral session
5. (todo) Populate foraging tables (especially model fitting)
'''

import subprocess
import sys
import logging
import json
import datetime
import os.path
import pandas as pd

#=============================   Change me!!! ===============================
# Address of remote training rig PCs
rigs = [
    #{'local': 'AIND-Tower-1', 'remote': R'\\10.128.37.23\Documents\Pybpod\Projects', 'user_name': 'labadmin', 'passcode': 'cupcake'},
    {'local': 'AIND-Tower-1', 'remote': R'\\10.128.37.23\Documents\Pybpod', 'user_name': 'labadmin', 'passcode': 'cupcake'},
    {'local': 'AIND-Tower-2', 'remote': R'\\10.128.200.146\Users\labadmin\Documents\foraging_projects\Projects', 'user_name': 'labadmin', 'passcode': 'cupcake'},
    {'local': 'AIND-Tower-2', 'remote': R'\\10.128.200.146\Users\labadmin\Documents\PyBpod', 'user_name': 'labadmin', 'passcode': 'cupcake'},
    {'local': 'AIND-Tower-3', 'remote': R'\\10.128.203.121\Users\labadmin\Documents\PyBpod\Projects\DO_NOT_USE', 'user_name': 'labadmin', 'passcode': 'cupcake'}, #old loc
    {'local': 'AIND-Tower-3', 'remote': R'\\10.128.203.121\Users\labadmin\Documents\PyBpod\Projects', 'user_name': 'labadmin', 'passcode': 'cupcake'},
    {'local': 'AIND-Ephys-Han', 'remote': R'\\10.128.54.220\Users\Han2\Documents\Pybpod\Projects', 'user_name': 'Han2', 'passcode': 'cupcake'},
    {'local': 'AIND-Tower-4', 'remote': R'\\10.128.37.31\Users\aind_behavior\Documents\PyBpod', 'user_name': 'aind_behavior', 'passcode': 'TraiNINGlab587!'},
    {'local': 'AIND-Tower-5', 'remote': R'\\10.128.41.7\Users\aind_behavior\Documents\PyBpod', 'user_name': 'aind_behavior', 'passcode': 'TraiNINGlab587!'},   
    {'local': 'AIND-Tower-6', 'remote': R'\\10.128.37.30\Users\aind_behavior\Documents\PyBpod', 'user_name': 'aind_behavior', 'passcode': 'TraiNINGlab587!'},
    {'local': 'AIND-Tower-7', 'remote': R'\\10.128.41.8\Users\aind_behavior\Documents\PyBpod', 'user_name': 'aind_behavior', 'passcode': 'TraiNINGlab587!'},
]

# Solve connection bugs: https://stackoverflow.com/questions/24933661/multiple-connections-to-a-server-or-shared-resource-by-the-same-user-using-more

# Path of dj config
dj_root = R'D:\Han_Sync\Svoboda\Scripts\map-ephys'

# Root path to place bpod projects (should be the same as in dj_local_conf.json)
behavioral_root = R'F:\Data_for_ingestion\Foraging_behavior\Behavior_rigs'

isilon_root = R'\\allen\programs\aind\workgroups\ephys\HanHou\BehaviorRaw'

# In Teams - Foraging Behavior - Files - Foraging training - metadata - ... (more options) - "Add shortcut to OneDrive", and locate the folder on your local PC
local_sharepoint_sync_folder = R'C:/Users/admin/OneDrive - Allen Institute/metadata'
remote_meta_file_animal = R'Surgery, water restriction and training.xlsx'
remote_meta_file_lab = R'Lab metadata.xlsx'

# Log files
pipeline_log = R'F:\Data_for_ingestion\Foraging_behavior\Behavior_rigs\behavior_automation.txt'  # Simplified log for this code
copy_log = R'F:\Data_for_ingestion\Foraging_behavior\Behavior_rigs\robocopy.txt'  # Full log for robocopy
#==========================================================================

# Get targeted folder from datajoint
with open(dj_root + '\dj_local_conf.json') as f:
    dj_conf = json.load(f)
local_meta_dir_animal = dj_conf['custom']['behavior_bpod']['meta_dir']
local_meta_dir_lab = dj_conf['custom']['behavior_bpod']['meta_lab_dir']

logging.basicConfig(# filename=pipeline_log, 
                    level=logging.INFO,
                    format='%(asctime)s %(name)s %(levelname)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    #filemode='a',
                    handlers=[logging.StreamHandler(), logging.FileHandler(pipeline_log)])

log = logging.getLogger(__name__)

def sync_behavioral_folders():
    
    for rig in rigs:
        summary_start = False
        command = fR'''net use {rig['remote']} /u:{rig['user_name']} {rig['passcode']}&&'''\
                  fR'''robocopy  {rig['remote']} {behavioral_root}\{rig['local']} /e /xx /XD "experiments_exported" "sessions_bak" "old_stuff" /xj /xjd /mt /np /Z /W:1 /R:5 /tee /fft /log+:{copy_log}&&'''\
                  fR'''net use {rig['remote']} /d'''         
                   ##fR'''net use {rig['remote']} /u:{rig['user_name']} {rig['passcode']}&&'''\          
                
        log.info('')
        log.info(f'''===== Sync "{rig['local']}" from {rig['remote']} ===== ''')
        proc = subprocess.Popen(command,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT,
                                shell=True,
                                universal_newlines=True  # Output to console
                                )
        
        # Show progress in the console
        for line in proc.stdout:
            sys.stdout.write(line)  # All output to console
            if 'Total' in line:     # Only summary part to log file
                summary_start = True
                log.info('')
            if summary_start:
                log.info(line.rstrip('\n'))
       
        exitcode = proc.wait() # Essential for robocopy
        log.info(f'Done with exitcode = {exitcode}\n')
        
def backup_to_isilon():
    command = fR'''robocopy F:\Data_for_ingestion\Foraging_behavior {isilon_root} /e /XD "experiments_exported" /xj /mt /Z /W:1 /R:5 /tee /fft'''
    proc = subprocess.Popen(command,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        shell=True,
                        universal_newlines=True  # Output to console
                        )
    for line in proc.stdout:
        sys.stdout.write(line)  # All output to console

    proc.wait()
        
        
def remote_meta_to_csvs(remote_notebook_name, local_meta_dir, transposed = False):
    # Save remote notebook to local csv files, one for each sheet
    # Update anyway (gspread with v4 api is not easy to get the last modi`fied time)

    file_path = os.path.join(local_sharepoint_sync_folder, remote_notebook_name)
    excel_file = pd.ExcelFile(file_path)
    last_mod = datetime.datetime.fromtimestamp(os.path.getmtime(file_path))
    log.info(f'Last modified: {last_mod}')
    
    sheet_names = excel_file.sheet_names

    for sheet in sheet_names:
        try:
            df_wr = pd.read_excel(excel_file, sheet)
                    
            if type(df_wr) == pd.DataFrame:
                df_wr.to_csv(os.path.join(local_meta_dir, f'{sheet}.csv')) 
                log.info(f'  csv saved: {sheet} in {remote_notebook_name}')
            else:
                log.info(f'  csv failed to save: {sheet} in {remote_notebook_name}')
        except:
            log.info(f'  csv failed to save: {sheet} in {remote_notebook_name}')
            
            
def ingest_behavior():
    os.chdir(dj_root)
    from pipeline.shell import load_meta_foraging, ingest_foraging_behavior, logsetup
    logsetup('INFO')
    
    try:
        load_meta_foraging()
    except:
        pass
    
    ingest_foraging_behavior()
                
if __name__ == '__main__':
    
  
    # Copy SharePoint spreadsheets meta info to local
    log.info('Updating surgery and WR metadata from Sharepoint sync folder...')
    remote_meta_to_csvs(remote_meta_file_animal, local_meta_dir_animal)
    log.info('done!\n')
        
    log.info('Updating Lab metadata from Sharepoint sync folder...')
    remote_meta_to_csvs(remote_meta_file_lab, local_meta_dir_lab)
    log.info('done!\n')
    
    # Copy behavioral folders from remote PCs to local
    sync_behavioral_folders()
    
    # Ingest behavior to datajoint
    ingest_behavior()
    
    # Backup behavior to isilon
    # backup_to_isilon()
