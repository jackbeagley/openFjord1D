# Load the .cpython module in build
import sys
sys.path.append('./build')
import numpy as np

import basin_diffusion as bd
import matplotlib.pyplot as plt

import cmocean as cmo
import gsw

import pandas as pd
import itertools

# https://pypi.org/project/matplotlib-label-lines/
from labellines import labelLine, labelLines

def flip(items, ncol):
    return itertools.chain(*[items[i::ncol] for i in range(ncol)])

# n_z = 8
# n_t = 2
n_z=1000
n_t=1000

N_DIFFUSIVITIES = 5

sa = np.empty((n_z, n_t), order='F')
tis = np.empty_like(sa)
rho = np.full_like(sa, np.nan)
kappa = np.full_like(sa, np.nan)
derived = np.empty((n_z, n_t, bd.indices.n_derived), order='F')

# z = np.linspace(-130, -430, n_z) # Malaspina Reach
z = np.linspace(-55, -127, n_z) # Deep Cove
# z = np.linspace(0.0, -10.0, n_z)
# sa[:, 0] = np.linspace(35.32, 35.24, n_z)
# tis[:, 0] = np.linspace(13.0, 12.2, n_z)
sa[:, 0] = np.linspace(35.32, 35.24, n_z)
tis[:, 0] = np.linspace(13.0, 12.2, n_z)

# Deep Cove densest
sa[:, 0] = np.linspace(35.33, 35.33, n_z)
tis[:, 0] = np.linspace(12.8, 12.8, n_z)

# hypso_df = pd.read_csv('malaspina_hypsography.csv')
hypso_df = pd.read_csv('deep_cove_south_hypsography.csv')

# resample the hypsography by interpolating at specified z values
a = np.interp(z, hypso_df['z'], hypso_df['a'])

# a = np.linspace(30000*1500, 10000*500, n_z)
# a = np.linspace(10.0, 10.0, n_z)

year_seconds = 60*60*24*365
month_start = 0
month_end = 12

t_start = month_start * year_seconds / 12
t_end = month_end * year_seconds / 12

t = np.linspace(t_start, t_end, n_t)

dt = (t_end - t_start) / (n_t - 1)
dz = (z[-1] - z[0]) / (n_z - 1)

# approximated from Moana
# sa_bc = sa_offshore_f(t/year_seconds*365)
# tis_bc = tis_offshore_f(t/year_seconds*365)
# sa_bc = np.full_like(t, 35.23)
# tis_bc = np.full_like(t, 12.9)

# lightest from intermediate water Deep Cove
sa_bc = np.full_like(t, 34.39)
tis_bc = np.full_like(t, 15.7)


diffusivities = np.logspace(-5, -3, N_DIFFUSIVITIES)

# Figure for heatmaps
fig = plt.figure(figsize=(8, 8))
fig_poster = plt.figure(figsize=(22.5/2.54, 9/2.54))

plt.rcParams.update({'font.size': 16})


# create a gridspec
gs_heatmap = fig.add_gridspec(N_DIFFUSIVITIES, 2)

# Plot monthly density profiles
fig_heatmap = plt.figure(figsize=(8, 8))

# create a gridspec
gs = fig_heatmap.add_gridspec(2, 3)
gs_poster = fig_poster.add_gridspec(1, 4)

ax_sigma0 = fig.add_subplot(gs[0, 0])

ax_sigma0.set_xlabel('Sigma-$\\theta$ Density \n/ $kgm^{-3}$')
ax_sigma0.set_ylabel('Depth / m')

ax_sal = fig.add_subplot(gs[0, 1], sharey=ax_sigma0)

ax_sal.set_xlabel('Absolute Salinity \n/ gkg$^{-1}$')
# ax_sal.set_ylabel('Depth / m')

ax_temp = fig.add_subplot(gs[0, 2], sharey=ax_sigma0)

ax_temp.set_xlabel('Temperature \n/ $^{\circ}$C')
# ax_temp.set_ylabel('Depth / m')

ax_deviation = fig.add_subplot(gs[1, :])
ax_deviation = fig_poster.add_subplot(gs_poster[0, 0:3])
ax_deviation.set_xlabel('Time / months')
ax_deviation.set_ylabel('Likeness / %')
# ax_deviation.set_ylabel('Likeness to Intermediate\nWater / %')

ax_initial_final = fig_poster.add_subplot(gs_poster[0, 3])
ax_initial_final.set_xlabel('Temperature / $^{\circ}$C')
ax_initial_final.set_ylabel('Depth / $m$')
ax_initial_final.yaxis.set_label_position("right")
ax_initial_final.yaxis.tick_right()
ax_initial_final.invert_yaxis()

ax_sigma0.plot(gsw.sigma0(sa[:, 0], tis[:, 0]), z, color='k', label='Initial', zorder=100)
ax_sal.plot(sa[:, 0], z, color='k', label='Initial', zorder=100)
ax_temp.plot(tis[:, 0], z, color='k', label='Initial', zorder=100)

alphas = np.linspace(0.25, 1.0, len(diffusivities))[::-1]

final_display_index = 1

for i, kappa_const in enumerate(diffusivities):
    bd.basin.compute(z, a, sa, tis, derived, t, sa_bc, tis_bc, bd.computation_modes.normal, bd.diffusivity_modes.constant, [kappa_const, np.nan])

    # Calculate sigma0
    sigma0 = gsw.sigma0(sa, tis)

    max_sa_deviation = np.max(np.abs(sa[:, 0] - sa_bc[0]))
    max_tis_deviation = np.max(np.abs(tis[:, 0] - tis_bc[0]))
    
    sa_range = np.max(sa, axis=0) - np.min(sa, axis=0)
    tis_range = np.max(tis, axis=0) - np.min(tis, axis=0)

    # calculate the sa and tis deviations as a percentage of maximum
    sa_deviation_percent = (sa - sa_bc[0]) / max_sa_deviation * 100
    tis_deviation_percent = (tis - tis_bc[0]) / max_tis_deviation * 100

    sa_range_norm = sa_range / np.max(sa_range)
    tis_range_norm = tis_range / np.max(tis_range)

    t_months = t/year_seconds*12

    tt_months, zz = np.meshgrid(t_months, z)

    heatmap_extent = [np.min(t_months), np.max(t_months), np.min(z), np.max(z)]

    label = '$\kappa={}$'.format(np.format_float_scientific(kappa_const, precision=2, exp_digits=1))

    ax_deviation.plot(t_months, 100 - np.max(np.abs(tis_deviation_percent), axis=0), color='r', alpha=alphas[i], label=label)
    # ax_deviation.plot(t_months, sa_range_norm*100, color='r', linestyle='--', alpha=alphas[i])

    axis_deviation = 100 - np.max(np.abs(tis_deviation_percent), axis=0)

    # add text to the right of each curve with the deviation value
    if i != 3:
        ax_deviation.text(t_months[-1], axis_deviation[-1], '{:.1f}%'.format(axis_deviation[-1]) , horizontalalignment='center', verticalalignment='center', bbox=dict(facecolor='white', alpha=0.75, edgecolor='none'))

    if i == final_display_index:
        percentage_tis_map_points = [
            [tis[0, 0], 100.0],
            [np.max(tis), 0.0]
        ]

        # percentage to temperature
        f = lambda x: (x - (percentage_tis_map_points[0][1])) / (percentage_tis_map_points[1][1] - percentage_tis_map_points[0][1]) * (percentage_tis_map_points[1][0] - percentage_tis_map_points[0][0]) + percentage_tis_map_points[0][0]
        # temperature to percentage
        f_inv = lambda x: (x - percentage_tis_map_points[0][0]) / (percentage_tis_map_points[1][0] - percentage_tis_map_points[0][0]) * (percentage_tis_map_points[1][1] - percentage_tis_map_points[0][1]) + percentage_tis_map_points[0][1]

        ax_deviation_2 = ax_initial_final.secondary_xaxis('top', functions=(f_inv, f))
        ax_deviation_2.set_xlabel('Likeness / %')

        ax_initial_final.plot(tis[:, 0], -z, color='k')#, label='Initial')
        ax_initial_final.plot(tis[:, -1], -z, color='r', label=(label))

        ax_initial_final.text(np.mean(tis[:, 0]), -np.mean(z), 'Initial', fontsize=14, rotation=90, horizontalalignment='center', verticalalignment='center', bbox=dict(facecolor='white', alpha=0.75, edgecolor='none'))
        
        labelLines(ax_initial_final.get_lines(), align=True, fontsize=14, xvals=(14.7, 14.7))

        # ax_initial_final.legend(markerscale=0.75)

    ax_sigma0.plot(sigma0[:, -2], z, color='r', alpha=alphas[i])

    ax_sal.plot(sa[:, -2], z, color='r', alpha=alphas[i], label=label)

    ax_temp.plot(tis[:, -2], z, color='r', alpha=alphas[i], label=label)

    # plot the heatmaps
    ax_sa_heatmap = fig_heatmap.add_subplot(gs_heatmap[i, 0])
    ax_temp_heatmap = fig_heatmap.add_subplot(gs_heatmap[i, 1])

    ax_sa_heatmap.imshow(sa[::-1, :], extent=heatmap_extent, aspect='auto', origin='lower', cmap=cmo.cm.haline)
    ax_temp_heatmap.imshow(tis[::-1, :], extent=heatmap_extent, aspect='auto', origin='lower', cmap=cmo.cm.thermal)

    ax_sa_heatmap.figure.colorbar(ax_sa_heatmap.images[0], ax=ax_sa_heatmap, label='Absolute Salinity / gkg$^{-1}$')
    ax_temp_heatmap.figure.colorbar(ax_temp_heatmap.images[0], ax=ax_temp_heatmap, label='Temperature / $^{\circ}$C')

handles, labels = ax_sal.get_legend_handles_labels()
leg = ax_sal.legend(flip(handles, 3), flip(labels, 3), loc='lower center', bbox_to_anchor=(0.5, 1.1), ncols=3)

# leg = .legend()
leg.set_in_layout(False)

for ax in [ax_sigma0, ax_sal, ax_temp]:
    ax.grid('on')

for ax in [ax_sal, ax_temp]:
    plt.setp(ax.get_yticklabels(), visible=False)

labelLines(ax_deviation.get_lines(), align=True, fontsize=14, xvals=(11, 1))

gs.tight_layout(fig, rect=[0, 0, 1, 0.9])
gs_heatmap.tight_layout(fig_heatmap, rect=[0, 0, 1, 0.9])

## make heatmap of evolution

fig_poster.savefig('diffusivity_comparison.pdf', bbox_inches='tight')

plt.show()