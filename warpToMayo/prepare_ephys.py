# -*- coding: utf-8 -*-
"""
Created on Sun Mar 14 21:05:56 2021

@author: houha
"""

from pathlib import Path
from atlaselectrophysiology.extract_files import extract_data

# Path to KS2 output
ks_path = Path(r'F:\catGT\HH102\catgt_HH102S07_C03P01_g0\HH102S07_C03P01_g0_imec0\imec0_ks2_orig')  # Using the orginal ks2 output is recommended!

# Path to raw ephys data
ephys_path = Path(r'F:\catGT\HH102\catgt_HH102S07_C03P01_g0\HH102S07_C03P01_g0_imec0')

# Save path
out_path = Path(r'F:\catGT\HH102\catgt_HH102S07_C03P01_g0\HH102S07_C03P01_g0_imec0\alf')

extract_data(ks_path, ephys_path, out_path)