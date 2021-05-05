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
styles_single_color = [['-', 2], ['-', 1], ['--', 1], [':', 1], ['-.', 0.5]]
styles_multiple_colors = [[[1.0000, 0.0078, 0.0078], '-', 2],
                          [[0.03137, 0.5725, 0.9882], '-', 2],
                          [[1.0000, 0.0078, 0.0078], '--', 1],
                          [[0.03137, 0.5725, 0.9882], '--', 1],
                          [[0.87, 0.87, 0.87], '--', 1],
                          ]


def _make_default_colormap():
    """Return the default colormap, with custom first colors."""
    colormap = np.array(cc.glasbey_bw_minc_20_minl_30)
    # Reorder first colors.
    colormap[[0, 1, 2, 3, 4, 5]] = colormap[[3, 0, 4, 5, 2, 1]]
    # Replace first two colors.
    colormap[0] = [0.03137, 0.5725, 0.9882]
    colormap[1] = [1.0000, 0.0078, 0.0078]
    return colormap


def load_bitcode_mat():  # Directly load and parse mat file

    # Load bitcode.mat
    mat = scipy.io.loadmat(glob.glob('*bitcode.mat')[0])
    dig_marker_per_trial = mat['digMarkerPerTrial']
    # Trial types
    ignore_trials = np.all(np.isnan(dig_marker_per_trial[:, [CHOICEL_, CHOICER_]]), 1)
    reward_trials = ~np.isnan(dig_marker_per_trial[:, REWARD_])
    noreward_trials = ~ignore_trials & ~reward_trials
    choiceL_trials = ~np.isnan(dig_marker_per_trial[:, CHOICEL_])
    choiceR_trials = ~np.isnan(dig_marker_per_trial[:, CHOICER_])
    choice_times = np.nanmean(dig_marker_per_trial[:, [CHOICEL_, CHOICER_]], 1)
    events = {}

    # ----- Define events times ------
    # 1. Choice_direction
    # events['choice_direction'] = {'Left': mat['choiceL'], 'Right': mat['choiceR']}

    # 2. Choice_outcome
    # choice_reward = choice_times[reward_trials]
    # choice_noreward = choice_times[noreward_trials]
    # events['choice_outcome'] = {'reward': choice_reward, 'no_reward': choice_noreward}

    # 1. Gocue_outcome
    # gocue_reward = dig_marker_per_trial[reward_trials, GOCUE_]
    # gocue_noreward = dig_marker_per_trial[noreward_trials, GOCUE_]
    gocue_L_reward = dig_marker_per_trial[reward_trials & choiceL_trials, GOCUE_]
    gocue_R_reward = dig_marker_per_trial[reward_trials & choiceR_trials, GOCUE_]
    gocue_L_noreward = dig_marker_per_trial[noreward_trials & choiceL_trials, GOCUE_]
    gocue_R_noreward = dig_marker_per_trial[noreward_trials & choiceR_trials, GOCUE_]
    gocue_ignore = dig_marker_per_trial[ignore_trials, GOCUE_]
    # events['gocue_direction_outcome'] = {'reward': gocue_reward, 'no_reward': gocue_noreward, 'ignore': gocue_ignore}
    events['gocue_direction_outcome'] = {'L_reward': gocue_L_reward, 'R_reward': gocue_R_reward,
                                         'L_noreward': gocue_L_noreward, 'R_noreward': gocue_R_noreward,
                                         'ignore': gocue_ignore}
    # 2. Choice_direction_outcome
    choice_L_reward = choice_times[reward_trials & choiceL_trials]
    choice_R_reward = choice_times[reward_trials & choiceR_trials]
    choice_L_noreward = choice_times[noreward_trials & choiceL_trials]
    choice_R_noreward = choice_times[noreward_trials & choiceR_trials]
    events['choice_direction_outcome'] = {'L_reward': choice_L_reward, 'R_reward': choice_R_reward,
                                          'L_noreward': choice_L_noreward, 'R_noreward': choice_R_noreward}

    # iti_reward = dig_marker_per_trial[reward_trials, ITI_]
    # iti_noreward = dig_marker_per_trial[noreward_trials, ITI_]
    # iti_ignore = dig_marker_per_trial[ignore_trials, ITI_]
    # events['iti_outcome'] = {'reward': iti_reward, 'no_reward': iti_noreward, 'ignore': iti_ignore}

    # 3. ITI_choice*outcome
    iti_L_reward = dig_marker_per_trial[reward_trials & choiceL_trials, ITI_]
    iti_R_reward = dig_marker_per_trial[reward_trials & choiceR_trials, ITI_]
    iti_L_noreward = dig_marker_per_trial[noreward_trials & choiceL_trials, ITI_]
    iti_R_noreward = dig_marker_per_trial[noreward_trials & choiceR_trials, ITI_]
    events['iti_direction_outcome'] = {'L_reward': iti_L_reward, 'R_reward': iti_R_reward,
                                       'L_noreward': iti_L_noreward, 'R_noreward': iti_R_noreward}

    return events


class PSTHView(ManualClusteringView):
    plot_canvas_class = PlotCanvasMpl  # use matplotlib instead of OpenGL (the default)

    def __init__(self, c=None):
        """features is a function (cluster_id => Bunch(data, ...)) where data is a 3D array."""
        super(PSTHView, self).__init__()
        self.controller = c
        self.model = c.model
        self.supervisor = c.supervisor
        self.cmap = _make_default_colormap()
        self.nodata = False

        try:
            self.events_all = load_bitcode_mat()
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
            if ((d - 1) // self.n_event_plots + 1) % 2 == 0:  # This is very awkward...
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

        for i, d in enumerate(cluster_ids):  # For each cluster
            psth_max_ylim = 0

            for j, (event_plot_name, event_plot_events) in enumerate(self.events_all.items()):  # For each event plot
                # Prepare
                raster_top = 0
                axis_psth = axis_list[i * (2 * self.n_event_plots) + j]
                axis_raster = axis_list[i * (2 * self.n_event_plots) + j + self.n_event_plots]

                axis_psth.axvline(x=0, color='white', alpha=.5)

                for k, (events_name, events_this) in enumerate(event_plot_events.items()):  # For each line in each plot
                    # Get data
                    rasters, activity, yrast, ntrials, nevents = self.get_spikes(d, events_this)
                    hist, bins = np.histogram(activity, weights=np.ones(nevents) * (50 / ntrials), range=(-5, 5),
                                              bins=250)

                    if len(cluster_ids) > 1:  # If more than one clusters are selected, only use the cluster's color
                        this_style, this_lw = styles_single_color[k % len(styles_single_color)]
                        this_color = self.cmap[i]
                    else:  # If only one cluster is selected, use Svoboda lab's convention (e.g., red for L, blue for R)
                        this_color, this_style, this_lw = styles_multiple_colors[k % len(styles_multiple_colors)]

                    psth_style = {'color': this_color, 'linestyle': this_style, 'lw': this_lw}
                    raster_style = {'color': [this_color]}

                    # Plot psth
                    axis_psth.plot(bins[:-1], hist, label=events_name, **psth_style)
                    psth_max_ylim = max(psth_max_ylim, max(hist) * 1.2)
                    # axis_psth.set_xticks(np.linspace(-2, 3, 9))

                    # Plot raster
                    axis_raster.scatter(rasters, yrast + raster_top,
                                        marker='o', s=np.ones(len(rasters)) * .4, alpha=0.8, **raster_style)
                    raster_top += ntrials
                    axis_raster.axhline(y=raster_top, color='white', alpha=.5)

                # Axis settings
                axis_psth.set_xlim(left=-2, right=3)
                axis_psth.set_title((f'#{d}: ' if j == 0 else '') + event_plot_name, fontsize=20)
                axis_psth.legend(loc='upper left')
                # if j > 0: axis_psth.set_yticks([])

                # psth_max_ylim = max(psth_max_ylim, max(axis_psth.get_ylim()))

                axis_raster.axvline(x=0, color='white', alpha=.5, lw=1)
                axis_raster.set_ylim(bottom=0, top=raster_top)
                axis_raster.invert_yaxis()
                axis_raster.set_xlabel(f'Time from {event_plot_name.split("_")[0]} (sec)', fontsize=15)
                self.fix_axis(axis_psth, 15)
                self.fix_axis(axis_raster, 15)

            # Same y axes for psth plots
            for j, _ in enumerate(self.events_all.items()):
                axis_psth = axis_list[i * (2 * self.n_event_plots) + j]
                axis_psth.set_ylim(bottom=0, top=psth_max_ylim)

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


class PSTHPlugin(IPlugin):
    def attach_to_controller(self, controller):
        def create_event_view():
            """A function that creates and returns a view."""
            return PSTHView(c=controller)

        controller.view_creator['PSTHView'] = create_event_view
