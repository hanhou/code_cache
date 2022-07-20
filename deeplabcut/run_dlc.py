#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 10 16:04:37 2019
@author: alex, han
"""

import os
from pathlib import Path
import re
import deeplabcut
import random
import multiprocessing as mp

from matplotlib.animation import FuncAnimation


def getsubfolders(folder):
    ''' returns list of subfolders '''
    return [os.path.join(folder,p) for p in os.listdir(folder) if os.path.isdir(os.path.join(folder,p))]


# ===================================================    
if_analyze = False
to_generate_labeled_video = 10  # Number of example labeled videos to generate in each folder

prefix=R'G:\DLC'
project='Foraging_Bot-Han_Lucas-2022-04-27'
shuffle=3

basepath='Z:\ephys\HanHou\Video' #data'
folder_must_have_pattern = R'.*_S.*'
folder_to_exclude = ['HH08_S01_20210812', 'HH08_S02_20210813', 'HH08_S03_20210814', 'HH08_S04_20210815']
file_must_include = 'bottom_face'

# ===================================================    

projectpath=os.path.join(prefix,project)
config=os.path.join(projectpath,'config.yaml')
subfolders=getsubfolders(basepath)


def dlc_batch_analyze():
    if if_analyze:
        for subfolder in subfolders:
            # if 'HH08' not in subfolder:
            #     continue
            
            print("Starting analyze data in:", subfolder)
            subsubfolders=getsubfolders(subfolder)
            for subsubfolder in subsubfolders:
                if re.search(folder_must_have_pattern, subsubfolder) is None: continue
                if subsubfolder in folder_to_exclude: continue
                
                print("Starting analyze data in:", subsubfolder)
                for vtype in ['.avi']:
                    deeplabcut.analyze_videos(config, [subsubfolder], must_include=file_must_include, shuffle=shuffle,videotype=vtype,save_as_csv=True, allow_growth=True)


def dlc_batch_generate_labeled_videos():
    to_generate_label = []
    
    # Get all videos to process
    for subfolder in subfolders:
        subsubfolders=getsubfolders(subfolder)
        for subsubfolder in subsubfolders:
            if re.search(folder_must_have_pattern, subsubfolder) is None: continue
            if subsubfolder in folder_to_exclude: continue

            all_videos = [str(v) for v in Path(subsubfolder).glob('*.avi') if file_must_include in str(v)]
            to_generate_label.extend(random.sample(all_videos, min(len(all_videos), to_generate_labeled_video)))
            
    # Do parallel processing            
    deeplabcut.create_labeled_video(config, shuffle=shuffle, 
                                    videos=to_generate_label, 
                                    save_frames=True, # High resolution 
                                    pool=pool)
        
            
if __name__ == '__main__':    
    if if_analyze:
        dlc_batch_analyze()
    
    if to_generate_labeled_video > 0:
        cores = int(mp.cpu_count()) - 1  # Auto core number selection
        pool = mp.Pool(processes=cores)  # For Windows, must be put inside "if __name__ == '__main__'" instead of any subfunctions
        
        dlc_batch_generate_labeled_videos()

        if pool != '':
            pool.close()
            pool.join()