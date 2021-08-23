"""
Example code for plotting spike snippets exported from raw spikeGLX .bin data.

Spike waveforms were retrieve from raw .bin file (before any filtering in catGT etc.), using 
spike times provided by kilosort results. See here for details:
https://djoshea.github.io/neuropixel-utils/waveforms/#extracting-waveforms-via-kilosortdataset

The example spikes.npz file was from a Neuropixels 1.0 probe.
"""

#%%
import numpy as np
import matplotlib.pyplot as plt

# -- Get data --
spikes = np.load('spikes.npz')
time_ms = spikes['time_ms']  # [1, n_time_bin]
waveforms = spikes['waveform']  # [n_time_bin, n_spikes_per_cell (100), n_cells (10)]
n_cells = waveforms.shape[2]

# -- Plotting --
m = np.floor(np.sqrt(n_cells)).astype(int)
n = np.ceil(n_cells / m).astype(int)
fig, axs = plt.subplots(m, n, figsize=(m * 5, n * 3))
for i in range(n_cells):
    waveform = waveforms[:, :, i]
    ax = axs.flatten()[i]
    ax.plot(time_ms, waveforms[:, :, i], 'k', alpha=0.1)
    ax.set(title=f'unit #{i}')
    
    if i == (m - 1) * n:
        ax.set(xlabel='Time (ms)', ylabel='Raw trace (uV)')
    
