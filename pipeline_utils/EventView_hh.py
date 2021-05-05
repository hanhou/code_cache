# import from plugins/matplotlib_view.py
"""Show how to create a custom matplotlib view in the GUI."""

from phy import IPlugin
from phy.cluster.views import ManualClusteringView  # Base class for phy views
from phy.plot.plot import PlotCanvasMpl  # matplotlib canvas
from numpy import genfromtxt
import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import colorcet as cc
import numpy as np
import warnings
import time
import sys
import os
import scipy.io
import glob

warnings.filterwarnings("ignore")

axis_list = []
plot_handles = []

STRIG_, GOCUE_, CHOICEL_, CHOICER_, REWARD_, ITI_ = 0, 1, 2, 3, 4, 5
line_styles = ['-', '--', ':', '-.']


def _make_default_colormap():
    """Return the default colormap, with custom first colors."""
    colormap = np.array(cc.glasbey_bw_minc_20_minl_30)
    # Reorder first colors.
    colormap[[0, 1, 2, 3, 4, 5]] = colormap[[3, 0, 4, 5, 2, 1]]
    # Replace first two colors.
    colormap[0] = [0.03137, 0.5725, 0.9882]
    colormap[1] = [1.0000, 0.0078, 0.0078]
    return colormap


class EventView(ManualClusteringView):
    plot_canvas_class = PlotCanvasMpl  # use matplotlib instead of OpenGL (the default)

    def __init__(self, c=None):
        """features is a function (cluster_id => Bunch(data, ...)) where data is a 3D array."""
        super(EventView, self).__init__()
        self.controller = c
        self.model = c.model
        self.supervisor = c.supervisor
        self.cmap = _make_default_colormap()
        self.nodata = False

        try:
            self.events_all = self.load_bitcode_mat()
            self.n_event_plots = len(self.events_all)
            # events = genfromtxt('events.csv', delimiter=',')
        except OSError:
            sys.stderr.write("EventView: bitcode_mat not found in: " + str(os.getcwd()) + "\n")
            # sys.stderr.write("EventView: events.csv not found in: " + str(os.getcwd()) + "\n")
            self.nodata = True

    def on_request_similar_clusters(self, cid=None):
        self.on_select()

    def on_select(self, cluster_ids=(), **kwargs):
        if self.nodata:
            return

        global axis_list, plot_handles
        cluster_ids = self.supervisor.selected
        self.cluster_ids = cluster_ids
        self.nclusts = len(cluster_ids)

        if axis_list:
            axis_diff = (self.nclusts - len(axis_list) // 2) * 2
            if axis_diff < 0:
                axis_list = axis_list[0:(len(axis_list)) + axis_diff]
                plot_handles = plot_handles[0:len(plot_handles) + axis_diff // 2]

        # We don't display anything if no clusters are selected.
        if not cluster_ids:
            return

        for i, d in enumerate(np.arange(start=1, stop=self.nclusts * 2 * self.n_event_plots + 1)):
            if ((d - 1) // self.n_event_plots + 1) % 2 == 0:   # This is very awkward...
                setattr(self, 'canvas.ax' + str(d), plt.subplot(2 * self.nclusts, self.n_event_plots, d,
                                                                sharex=axis_list[i - self.n_event_plots]))
            else:
                setattr(self, 'canvas.ax' + str(d), plt.subplot(2 * self.nclusts, self.n_event_plots, d))

            if (len(axis_list) - 1) < i:
                axis_list.append(getattr(self, 'canvas.ax' + str(d)))
            else:
                axis_list[i] = (getattr(self, 'canvas.ax' + str(d)))
            axis_list[i].cla()

        t1 = time.time()
        ttime = 0

        # import pdb;
        # pdb.set_trace()

        t = (time.time())

        for i, d in enumerate(cluster_ids):    # For each cluster
            for j, (event_plot_name, event_plot_events) in enumerate(self.events_all.items()):      # For each event plot
                raster_top = 0

                for k, (events_name, events_this) in enumerate(event_plot_events.items()):   # For each line in each plot
                    # Get data
                    rasters, activity, yrast, ntrials, nevents = self.get_spikes(d, events_this)
                    hist, bins = np.histogram(activity, weights=np.ones(nevents) * (50 / ntrials), range=(-5, 5), bins=250)

                    # Prepare
                    this_style = line_styles[k % len(line_styles)]
                    axis_psth = axis_list[i * (2 * self.n_event_plots) + j]
                    axis_raster = axis_list[i * (2 * self.n_event_plots) + j + self.n_event_plots]

                    # Plot psth
                    axis_psth.plot(bins[:-1], hist, color=self.cmap[i], linestyle=this_style)
                    axis_psth.axvline(x=0, color='white', alpha=.5)
                    axis_psth.set_xticks(np.linspace(-2, 3, 9))
                    axis_psth.set_xlim(left=-2, right=3)

                    # Plot raster
                    axis_raster.scatter(rasters, yrast + raster_top, c=[self.cmap[i]],
                                        marker='|', s=np.ones(len(rasters)) * .4, alpha=0.8)
                    axis_raster.axvline(x=0, color='white', alpha=.5)
                    raster_top += ntrials
                    axis_raster.axhline(y=raster_top, color='white', alpha=.5)

                    # Axis settings
                    axis_raster.set_ylim(bottom=0, top=raster_top)
                    axis_psth.set_ylim(bottom=0, top=None)
                    self.fix_axis(axis_psth, 10)
                    self.fix_axis(axis_raster, 10)


        print(ttime, 'plotting')
        t = time.time()
        # Use this to update the matplotlib figure.
        self.canvas.show()
        print(time.time() - t, 'update')
        print(time.time() - t1, 'total')
        return

    def fix_axis(self, ax, textsize):
        ax.tick_params(axis='x', labelsize=textsize)
        ax.tick_params(axis='y', labelsize=textsize)
        ax.xaxis.label.set_color('white')
        ax.tick_params(axis='x', colors='white')
        ax.tick_params(axis='y', colors='white')
        ax.grid(False)

    def get_spikes(self, clust, events):
        spikes = self.model.get_cluster_spikes(clust)
        spike_times = np.array(self.model.spike_times[spikes])
        rasters = np.array([])
        yrast = np.array([])
        activity = []
        last_ind = 0
        for i, d in enumerate(events):
            if d < np.amax(spike_times) and d > np.amin(spike_times):
                st = spike_times[last_ind:]
                temp1 = st - (d + 5)
                ind1 = np.abs(temp1).argmin()
                if temp1[ind1] > 0:
                    ind1 -= 1
                temp2 = st - (d - 5)
                ind2 = np.abs(temp2).argmin()
                if temp2[ind2] < 0:
                    ind2 += 1
                temp = st[ind2:ind1] - d
                last_ind = ind1
                rasters = np.append(rasters, temp)
                yrast = np.append(yrast, np.ones(len(temp)) * i)
                activity.extend(temp)
        return rasters, np.array(activity), yrast, np.amax(yrast), len(activity)

    def load_bitcode_mat(self):  # Directly load and parse mat file

        # Load bitcode.mat
        mat = scipy.io.loadmat(glob.glob('*bitcode.mat')[0])
        dig_marker_per_trial = mat['digMarkerPerTrial']
        # Trial types
        ignore_trials = np.all(np.isnan(dig_marker_per_trial[:, [CHOICEL_, CHOICER_]]), 1)
        reward_trials = ~np.isnan(dig_marker_per_trial[:, REWARD_])
        noreward_trials = ~ignore_trials & ~reward_trials
        events = {}

        # Define events times
        choice_L, choice_R = mat['choiceL'], mat['choiceR']
        events['choice_direction'] = {'Left': choice_L, 'Right': choice_R}

        gocue_reward = dig_marker_per_trial[reward_trials, GOCUE_]
        gocue_noreward = dig_marker_per_trial[noreward_trials, GOCUE_]
        gocue_ignore = dig_marker_per_trial[ignore_trials, GOCUE_]
        events['gocue_outcome'] = {'reward': gocue_reward, 'no_reward': gocue_noreward, 'ignore': gocue_ignore}

        choice_times = np.nanmean(dig_marker_per_trial[:, [CHOICEL_, CHOICER_]], 1)
        choice_reward = choice_times[reward_trials]
        choice_noreward = choice_times[noreward_trials]
        events['choice_outcome'] = {'reward': choice_reward, 'no_reward': choice_noreward}

        iti_reward = dig_marker_per_trial[reward_trials, ITI_]
        iti_noreward = dig_marker_per_trial[noreward_trials, ITI_]
        iti_ignore = dig_marker_per_trial[ignore_trials, ITI_]
        events['iti_outcome'] = {'reward': iti_reward, 'no_reward': iti_noreward, 'ignore': iti_ignore}

        return events


class EventPlugin(IPlugin):
    def attach_to_controller(self, controller):
        def create_event_view():
            """A function that creates and returns a view."""
            return EventView(c=controller)

        controller.view_creator['EventView'] = create_event_view
