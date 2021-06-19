# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 21:05:56 2021

@author: houha
"""

from pathlib import Path
from atlaselectrophysiology.extract_files import extract_data
import os


# root_path = r'i:\HH102\catgt_HH102S10_C04P02_g0\HH102S10_C04P02_g0_imec0'
# root_path = r'V:\Ingested\SC045\catgt_SC045_120920_g0\SC045_120920_g0_imec0'
root_path = r'v:\Ingested\SC045\catgt_SC045_120920_g0\SC045_120920_g0_imec2'
exclude_list = []

for root, subdirs, files in os.walk(root_path):
    for subdir in subdirs:
        if 'ks2_orig' in subdir:
            this_path = PATH_JOIN = os.path.join(root, subdir)
            if not any(ee in this_path for ee in exclude_list):
                print('Processing:', this_path, '...')

                # Path to KS2 output
                ks_path = Path(this_path)  # Using the original ks2 output is recommended!
                ephys_path = ks_path.parent  # Path to raw ephys data
                out_path = ephys_path.joinpath('alf')  # Save path

                if not out_path.exists():
                    out_path.mkdir()

                # extract_data(ks_path, ephys_path, out_path, if_ap_rmsmap=False)
                extract_data(ks_path, ephys_path, out_path)

