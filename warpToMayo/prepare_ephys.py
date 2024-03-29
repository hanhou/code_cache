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
root_paths =[R'G:\catGT\HH13', R'I:\HH13', R'J:', R'K:']
exclude_list = []

for root_path in root_paths:
    for root, subdirs, files in os.walk(root_path):
        for subdir in subdirs:
            if 'ks2' in subdir and 'orig' not in subdir:
                this_path = os.path.join(root, subdir)
                if not any(ee in this_path for ee in exclude_list):
                    print('\n\nProcessing:', this_path, '...')

                    # Path to KS2 output
                    ks_path = Path(this_path)  # Using the original ks2 output is recommended!
                    ephys_path = ks_path.parent  # Path to raw ephys data
                    out_path = ephys_path.joinpath('alf')  # Save path

                    if not out_path.exists():
                        out_path.mkdir()

                    extract_data(ks_path, ephys_path, out_path, max_length_in_sec=300)

