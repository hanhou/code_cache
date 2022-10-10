'''
1. From a root catGT folder, locate channel_location.json and xyz_pick from IBL GUI in each alf folder
--removed-- 2. Copy and rename them to the target folder so that they're ready for ingestion (format=2)
3. Flip LR for bug fix
'''

#%%
from pathlib import Path
import re, glob, shutil, json, os

#%%
root_catGT = R"z:\ephys\HanHou\catGT"
# root_for_ingestion = R"F:\Data_for_ingestion\MAP"

to_flip_LR = ['HH08', 'HH09']  # Temporary patch for not using the TrackFinder-flipped xyz_picks
to_adjust_minor_offset = ['HH13']  # Patch for 5700 is used in Dave's code but the actual LF flipping middle point should be 5739

#%%

all_alfs = glob.glob(root_catGT + R"\**\alf", recursive=True)

for alf in all_alfs:
    
    #%%
    
    # Get channel_locations file
    f_channel_location = list(Path(alf).glob("channel_locations*.json"))
    f_xyz_picks = list(Path(alf).glob("*xyz_picks*.json"))  
    if not (len(f_channel_location) and len(f_xyz_picks)):
        continue    
    
    # Get meta info
    meta = re.search(r'(.*)[\\|/](?P<h2o>[^_]*)_(?P<date>\w+)_(?P<g>g\d+)_(?P<imec>imec\d+)\\alf', alf)
    if meta is None:
        continue
    
    if meta['h2o'] not in (to_flip_LR + to_adjust_minor_offset):
        continue
    else:
        print(meta.group(0))
    
    # Fix channel_locations
    for f_ch in f_channel_location:    # Could be more than one shanks 
        with open(f_ch, 'r') as f:
            channel_locations = json.load(f)
        
        for pt in channel_locations.values():
            if 'x' in pt:
                if meta['h2o'] in to_flip_LR:
                    pt['x'] = - pt['x']
                elif meta['h2o'] in to_adjust_minor_offset:
                    pt['x'] = pt['x'] + 78  # In Dave's code, he used (2 * 5700 - Fiji), but it should be (2 * 5739 - Fiji)
        
        # Archive old file and save new file
        f_name = re.search(r'channel_(.*)', f_ch.name)
        os.rename(f_ch, f_ch.parent / ("channel_archive_" + f_name.groups()[0]))  # To break the pattern so that it won't be ingested
        with open(f_ch, 'w') as f:
            json.dump(channel_locations, f, indent=2, separators=(',', ': '))
            print(f'patched: {f_ch}')
        
    # Fix xyz picks (so that no mroe problems if re-run the GUI)
    for f_xyz in f_xyz_picks:    # Could be more than one shanks 
        with open(f_xyz, 'r') as f:
            xyz_picks = json.load(f)
        
        for xyz in xyz_picks['xyz_picks']:
            if meta['h2o'] in to_flip_LR:
                xyz[0] = - xyz[0]
            elif meta['h2o'] in to_adjust_minor_offset:
                xyz[0] = xyz[0] + 78  # In Dave's code, he used (2 * 5700 - Fiji), but it should be (2 * 5739 - Fiji)
        
        # Archive old file and save new file
        f_name = re.search(r'(.*)xyz_picks(.*)', f_xyz.name)
        os.rename(f_xyz, f_xyz.parent / (f_name.groups()[0] + "xyz_archive_picks" + f_name.groups()[1]))  # To break the pattern so that it won't be ingested
        with open(f_xyz, 'w') as f:
            json.dump(xyz_picks, f)
            print(f'patched: {f_xyz}')
        
       
            
        

    # # Find target folder
    # folder_for_ingestion = Path(root_for_ingestion) / meta['h2o'] / 'histology'
    # folder_for_ingestion.mkdir(parents=True, exist_ok=True)
    
    