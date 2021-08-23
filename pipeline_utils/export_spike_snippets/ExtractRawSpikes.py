#%%
import os
import scipy.io as spio
import numpy as np

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

spike_snippets = spio.loadmat('./spike_snippets.mat', squeeze_me=True)

np.savez('spikes', time_ms=spike_snippets['time_ms'], waveform=spike_snippets['spikes'])