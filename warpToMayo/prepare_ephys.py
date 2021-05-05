# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 21:05:56 2021

@author: houha
"""

from pathlib import Path
from atlaselectrophysiology.extract_files import extract_data
import os

PATH_JOIN = os.path.join(root, subdir)

root_path = r'v:\Ingested\SC061'
exclude_list = [r'SC045_120920_g0_imec0\imec0_ks2_orig']

for root, subdirs, files in os.walk(root_path):
    for subdir in subdirs:
        if 'ks2_orig' in subdir:
            this_path = PATH_JOIN
            if not any(ee in this_path for ee in exclude_list):
                print('Processing:', this_path, '...')

                # Path to KS2 output
                ks_path = Path(this_path)  # Using the original ks2 output is recommended!
                ephys_path = ks_path.parent  # Path to raw ephys data
                out_path = ephys_path.joinpath('alf')  # Save path

                extract_data(ks_path, ephys_path, out_path)
