"""
A lightweight util that calls CatGt commands recursively 
for all imec folders inside a user-defined root directory.
(so far tested for batch extracting lfp from both 1.0 and 2.0 probes)

Han @ Aug 2021
"""

import re
import os
import glob
import subprocess

# ========= Settings ==========
catgt_exe = r'D:\Han_Sync\Svoboda\Scripts\Ephys\CatGT_2.0\CatGT.exe'
data_root_path = r'h:' 

# If empty, output files will be stored in the original folders,
# otherwise, in conventional CatGt folders structure.
output_folder = r'f:\Test'

include_keywords = []  # If not empty, only search folders with the names, e.g., ['surface', '202101']
exclude_keywords = []  # If not empty, exclude folders with the names

# Shared CatGT command
max_secs = 100  # A minimal of 100 sec is recommended for better estimation of LFP surface
catgt_cmd = f'-lf -lfhipass=.1 -lflopass=300 -tshift -maxsecs={max_secs}'  
# =============================

# Get all folders with the pattern *g*_imec* recursively
imec_folders = glob.glob(data_root_path + '/**/*g*_imec*', recursive=True)
if output_folder: 
    os.makedirs(output_folder, exist_ok=True)

output_cmd = f'-out_prb_fld -dest={output_folder}' if output_folder else ''
n_probe = 0

for imec in imec_folders:
    if include_keywords and not any(inc in imec for inc in include_keywords):
        continue
    if any(ex in imec for ex in exclude_keywords):
        continue

    # Retrieve dir, run, g, prb
    m = re.search(r'(.*)[\\|/](.*)[\\|/](\w+)_g(\d+)_imec(\d+)', imec)
    if not m:
        print(f' !! Errors in parsing folder {imec}')
        continue
    dir, run, g, prb = m.group(1), m.group(3), m.group(4), m.group(5)

    # Retrieve t_max
    files = glob.glob(imec + '/*g*_t*.ap.bin') 
    if not files:    # Must have at least one .ap.bin file
        print(f' !! No ap.bin found for {imec}')
        continue
    ts = [int(re.search(r'\w+_g\d+_t(\d+).*', file).group(1)) for file in files]
    t_max = max(ts)

    # Do catGT
    cmd = f'{catgt_exe} -dir={dir} -prb_fld -run={run} -g={g} -prb={prb} -t=0,{t_max} {catgt_cmd} {output_cmd}'
    err = subprocess.call(cmd.split())  # Not using split(' ') to take care of more than one spaces between strings

    if not err:
        n_probe += 1
        print(f'{n_probe} done: {dir}, {run}, imec{prb}, g{g}, t_max = {t_max}')
    else:
        print(f' !! Errors in processing {imec}, t_max = {t_max}')

