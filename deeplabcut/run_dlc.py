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
if_analyze = True
to_generate_labeled_video_online = 10   # Number of example labeled videos to generate in each folder (during analyzing)
to_generate_labeled_video_offline = False  # Number of example labeled videos to generate in each folder (after analyzing)

# Project
prefix=R'G:\DLC'
project=R'Foraging_Bot-Han_Lucas-2022-04-27'
# project=R'SideView-Han_Lucas-2022-07-06'
shuffle=1

# To analyze
basepath='Z:\ephys\HanHou\Video' #data'
folder_must_have_pattern = R'.*_(S|photostim).*'
folder_to_exclude = ['HH08', 'HH09', ]
file_must_include = 'bottom_face'
#file_must_include = 'side_face'

# ===================================================    

projectpath=os.path.join(prefix,project)
config=os.path.join(projectpath,'config.yaml')
subfolders=getsubfolders(basepath)


def dlc_batch_analyze():
    if if_analyze:
        for subfolder in subfolders:
            subsubfolders=getsubfolders(subfolder)
            for subsubfolder in subsubfolders:
                if re.search(folder_must_have_pattern, subsubfolder) is None: 
                    # print(' no pattern matched,')
                    continue
                if any([s in subsubfolder for s in folder_to_exclude]):
                    # print(' excluded...')
                    continue
                
                # Skip if all videos have been analyzed before in a subsubfolder
                if len(list(Path(subsubfolder).glob(f'*{file_must_include}*.avi'))) == len(list(Path(subsubfolder).glob(f'*{file_must_include}*.csv'))):
                    print(f'       All {file_must_include} videos are done in {subsubfolder}, skipped!')
                    continue
                
                
                # Analyze videos
                print("\n\n========================Starting analyze data in: ", subsubfolder)
                for vtype in ['.avi']:
                    deeplabcut.analyze_videos(config, [subsubfolder], must_include=file_must_include, shuffle=shuffle,videotype=vtype,save_as_csv=True, allow_growth=True)
                    
                # Generate labeled video online
                try:
                    if to_generate_labeled_video_online:
                        # Do parallel processing        
                        print(f"\n-------------------\nGenerating {to_generate_labeled_video_online} labeled videos in:", subsubfolder)
    
                        to_generate_label = []
                        all_videos = [str(v) for v in Path(subsubfolder).glob('*.avi') if file_must_include in str(v)]
                        to_generate_label.extend(random.sample(all_videos, min(len(all_videos), to_generate_labeled_video_online)))

                        deeplabcut.create_labeled_video(config, shuffle=shuffle, 
                                                        videos=to_generate_label, 
                                                        save_frames=True, # High resolution 
                                                        pool=pool)
                except:
                    print('------------------- Error: something went wrong when generating labeled videos!! ----------------------')

                    

def dlc_batch_generate_labeled_videos_offline():
    to_generate_label = []
    
    # Get all videos to process
    for subfolder in subfolders:
        subsubfolders=getsubfolders(subfolder)
        for subsubfolder in subsubfolders:
            if re.search(folder_must_have_pattern, subsubfolder) is None: continue
            if subsubfolder in folder_to_exclude: continue

            all_videos = [str(v) for v in Path(subsubfolder).glob('*.avi') if file_must_include in str(v)]
            to_generate_label.extend(random.sample(all_videos, min(len(all_videos), to_generate_labeled_video_offline)))
            
    # Do parallel processing            
    deeplabcut.create_labeled_video(config, shuffle=shuffle, 
                                    videos=to_generate_label, 
                                    save_frames=True, # High resolution 
                                    pool=pool)
        
            
if __name__ == '__main__':    
    pool = ''
    
    if to_generate_labeled_video_online or to_generate_labeled_video_offline:
        cores = max(to_generate_labeled_video_online, to_generate_labeled_video_offline, 10)  # Auto core number selection
        pool = mp.Pool(processes=cores)  # For Windows, must be put inside "if __name__ == '__main__'" instead of any subfunctions
        
    if if_analyze:
        dlc_batch_analyze()
        
    if to_generate_labeled_video_offline:
        dlc_batch_generate_labeled_videos_offline()

    if pool != '':
        pool.close()
        pool.join()
