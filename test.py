# Load the .cpython module in build
import sys
sys.path.append('./build')
import numpy as np

import basin_diffusion as bd
import matplotlib.pyplot as plt

import cmocean as cmo
import gsw

import pandas as pd

# n_z = 8
# n_t = 2
n_z=1000
n_t=1000


sa_offshore_f = lambda d: 34.85 + 0.15 * np.cos((d - 120)/365 * 2*np.pi)
tis_offshore_f = lambda d: 12.5 + 2.0 * np.cos((d - 120)/365 * 2*np.pi)

sa = np.empty((n_z, n_t), order='F')
tis = np.empty_like(sa)
rho = np.full_like(sa, np.nan)
kappa = np.full_like(sa, np.nan)
derived = np.empty((n_z, n_t, bd.indices.n_derived), order='F')

z = np.linspace(-120, -430, n_z)
# z = np.linspace(0.0, -10.0, n_z)
sa[:, 0] = np.linspace(35.32, 35.24, n_z)
tis[:, 0] = np.linspace(13.0, 12.2, n_z)

hypso_df = pd.read_csv('malaspina_hypsography.csv')

# resample the hypsography by interpolating at specified z values
a = np.interp(z, hypso_df['z'], hypso_df['a'])

# a = np.linspace(30000*1500, 10000*500, n_z)
# a = np.linspace(10.0, 10.0, n_z)

year_seconds = 60*60*24*365
month_start = 4
month_end = 4 + 12

t_start = month_start * year_seconds / 12
t_end = month_end * year_seconds / 12

t = np.linspace(t_start, t_end, n_t)

dt = (t_end - t_start) / (n_t - 1)
dz = (z[-1] - z[0]) / (n_z - 1)

sa_bc = sa_offshore_f(t/year_seconds*365)
tis_bc = tis_offshore_f(t/year_seconds*365)

kappa_buoyancy_coefs = [2e-7, 1.6]
kappa_const = 1e-4

bd.basin.compute(z, a, sa, tis, derived, t, sa_bc, tis_bc, bd.computation_modes.closed_top, bd.diffusivity_modes.constant, [kappa_const, np.nan])

kappa = derived[:, :, bd.indices.kappa-1]
rho = derived[:, :, bd.indices.rho-1]
n_freq = derived[:, :, bd.indices.n_freq-1]

# Calculate sigma0
sigma0 = gsw.sigma0(sa, tis)

rho_grad = np.gradient(rho, axis=0) / dz

# create a mask for negative density gradients
instability_mask = rho_grad > 0.0
instability_mask_rgb = np.multiply(instability_mask[:, :, np.newaxis], [1.0, 0.0, 0.0])

# calculate Brunt-Vaisala frequency
N = np.sqrt(-(9.81 / rho) * rho_grad)

kappa_py = kappa_buoyancy_coefs[0] * N ** (-kappa_buoyancy_coefs[1])
kappa_py[np.isinf(kappa_py)] = np.nan

kappa_range = np.nanquantile(kappa, [0.01, 0.99])

t_months = t/year_seconds*12

tt_months, zz = np.meshgrid(t_months, z)

# Make a heatmap of salinity evolution
fig = plt.figure()
ax_sal = fig.add_subplot(3, 2, 1)
ax_sal.set_xlabel('Time / months')
ax_sal.set_ylabel('Depth / m')

heatmap_extent = [np.min(t_months), np.max(t_months), np.min(z), np.max(z)]

ax_sal.imshow(sa, cmap=cmo.cm.haline, aspect='auto', extent=heatmap_extent)

# Add a colorbar
cbar = ax_sal.figure.colorbar(ax_sal.images[0], ax=ax_sal, label='Absolute Salinity / gkg$^{-1}$')

# Add labelled contours
cont = ax_sal.contour(tt_months, zz, sa, levels=np.linspace(30, 40, 101), colors='k', linewidths=0.5)
ax_sal.clabel(cont, inline=True, fontsize=10, fmt='%1.1f')

# Make a heatmap of temperature evolution
ax_temp = fig.add_subplot(3, 2, 3)
ax_temp.set_xlabel('Time / months')
ax_temp.set_ylabel('Depth / m')

ax_temp.imshow(tis, cmap=cmo.cm.thermal, aspect='auto', extent=heatmap_extent)

# Add a colorbar
cbar = ax_temp.figure.colorbar(ax_temp.images[0], ax=ax_temp, label='Temperature / $^{\circ}$C')

# Add labelled contours
cont = ax_temp.contour(tt_months, zz, tis, levels=np.linspace(0, 20, 51), colors='k', linewidths=0.5)
ax_temp.clabel(cont, inline=True, fontsize=10, fmt='%1.1f')

# density
ax_sigma0 = fig.add_subplot(3, 2, 5)
ax_sigma0.set_xlabel('Time / months')
ax_sigma0.set_ylabel('Depth / m')  

ax_sigma0.imshow(sigma0, cmap=cmo.cm.dense, aspect='auto', extent=heatmap_extent)

# Add a colorbar
cbar = ax_sigma0.figure.colorbar(ax_sigma0.images[0], ax=ax_sigma0, label='Sigma0 Density / kg m$^{-3}$')

# show unstable regions as a transparent red
ax_sigma0.imshow(instability_mask_rgb, extent=heatmap_extent, aspect='auto', alpha=0.2)

# ax_density.plot(t_months, instability_depth, color='r', linewidth=1.0)

# Add labelled contours
cont = ax_sigma0.contour(tt_months, zz, sigma0, levels=np.linspace(20, 40, 201), colors='k', linewidths=0.5)
ax_sigma0.clabel(cont, inline=True, fontsize=10, fmt='%1.1f')

# diffusivity - kappa
ax_kappa = fig.add_subplot(3, 2, 2)
ax_kappa.set_xlabel('Time / months')
ax_kappa.set_ylabel('Depth / m')

ax_kappa.imshow(kappa, cmap=cmo.cm.amp, aspect='auto', extent=heatmap_extent, vmin=kappa_range[0], vmax=kappa_range[1])

# Add a colorbar
cbar = ax_kappa.figure.colorbar(ax_kappa.images[0], ax=ax_kappa, label='Diffusivity / m$^2$s$^{-1}$')

# Brunt-Vaisala frequency
ax_n = fig.add_subplot(3, 2, 4)
ax_n.set_xlabel('Time / months')
ax_n.set_ylabel('Depth / m')

ax_n.imshow(N, cmap=cmo.cm.speed, aspect='auto', extent=heatmap_extent)

cbar = ax_n.figure.colorbar(ax_n.images[0], ax=ax_n, label='Brunt-Vaisala frequency / s$^{-1}$')

# python derived diffusivity
ax_kappa_py = fig.add_subplot(3, 2, 6)
ax_kappa_py.set_xlabel('Time / months')
ax_kappa_py.set_ylabel('Depth / m')

ax_kappa_py.imshow(kappa_py, cmap=cmo.cm.amp, aspect='auto', extent=heatmap_extent, vmin=kappa_range[0], vmax=kappa_range[1])

cbar = ax_kappa_py.figure.colorbar(ax_kappa_py.images[0], ax=ax_kappa_py, label='Diffusivity / m$^2$s$^{-1}$')

plt.tight_layout()

# Plot monthly density profiles
fig = plt.figure()
ax_sigma0 = fig.add_subplot(1, 4, 1)
ax_sigma0.plot(sigma0[:, 0], z, color='k', label='Initial')
ax_sigma0.plot(sigma0[:, -2], z, color='r', label='Final')

ax_sigma0.set_xlabel('Sigma0 Density / kgm$^{-3}$')
ax_sigma0.set_ylabel('Depth / m')

ax_sigma0.legend()

ax_sal = fig.add_subplot(1, 4, 2, sharey=ax_sigma0)
ax_sal.plot(sa[:, 0], z, color='k', label='Initial')
ax_sal.plot(sa[:, -2], z, color='r', label='Final')

ax_sal.set_xlabel('Absolute Salinity / gkg$^{-1}$')
# ax_sal.set_ylabel('Depth / m')

ax_sal.legend()

ax_temp = fig.add_subplot(1, 4, 3, sharey=ax_sigma0)
ax_temp.plot(tis[:, 0], z, color='k', label='Initial')
ax_temp.plot(tis[:, -2], z, color='r', label='Final')

ax_temp.set_xlabel('Temperature / $^{\circ}$C')
# ax_temp.set_ylabel('Depth / m')

ax_temp.legend()

# plot the hypsographic curve
ax_hypso = fig.add_subplot(1, 4, 4, sharey=ax_sigma0)
ax_hypso.plot(a/1000000, z, color='k')

ax_hypso.set_xlabel('Area / km$^{2}$')
# ax_hypso.set_ylabel('Depth / m')

for ax in [ax_sigma0, ax_sal, ax_temp, ax_hypso]:
    ax.grid('on')

for ax in [ax_sal, ax_temp, ax_hypso]:
    plt.setp(ax.get_yticklabels(), visible=False)

plt.show()