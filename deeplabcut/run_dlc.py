#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 10 16:04:37 2019
@author: alex, han
"""

import os
import re
import deeplabcut

def getsubfolders(folder):
    ''' returns list of subfolders '''
    return [os.path.join(folder,p) for p in os.listdir(folder) if os.path.isdir(os.path.join(folder,p))]


# ===================================================    
prefix=R'G:\DLC'
project='Foraging_Bot-Han_Lucas-2022-04-27'
shuffle=3

basepath='Z:\ephys\HanHou\Video' #data'

# ===================================================    

projectpath=os.path.join(prefix,project)
config=os.path.join(projectpath,'config.yaml')


'''
Imagine that the data (here: videos of 3 different types) are in subfolders:
    /January/January29 ..
    /February/February1
    /February/February2
    
    etc.
'''

subfolders=getsubfolders(basepath)
for subfolder in subfolders: #this would be January, February etc. in the upper example
    # if 'HH08' not in subfolder:
    #     continue
    
    print("Starting analyze data in:", subfolder)
    subsubfolders=getsubfolders(subfolder)
    for subsubfolder in subsubfolders: #this would be Febuary1, etc. in the upper example...
        
        if re.search(R'.*_S.*', subsubfolder) is None:
            continue
        
        if subsubfolder in ['HH08_S01_20210812', 'HH08_S02_20210813', 'HH08_S03_20210814', 'HH08_S04_20210815']:
            continue
        
        print("Starting analyze data in:", subsubfolder)
        for vtype in ['.avi']:
            deeplabcut.analyze_videos(config,[subsubfolder],must_include='bottom_face', shuffle=shuffle,videotype=vtype,save_as_csv=True, allow_growth=True)