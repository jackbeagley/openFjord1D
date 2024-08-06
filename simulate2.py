# Load the .cpython module in build
import sys
sys.path.append('./build')
import numpy as np

import basin_diffusion as bd
import matplotlib.pyplot as plt

import cmocean as cmo
import gsw
cmo.cm
import pandas as pd

import xarray as xr

scenarios = {
	'Bran': {
		'ds': xr.open_dataset('offshore_doubtful_bran.nc'),
		},
	'Copernicus': {
		'ds': xr.open_dataset('offshore_doubtful_2015.nc'),
		},
	}

n_z = 1000

z = np.linspace(-55, -127, n_z)
hypso_df = pd.read_csv('deep_cove_north_hypsography.csv')
kappa_const = 1e-4  # Estimate for Deep Cove

#z = np.linspace(-130, -430, n_z)
#kappa_const = 1e-3  # Estimate for Malaspina Reach
#hypso_df = pd.read_csv('malaspina_hypsography.csv')

# resample the hypsography by interpolating at specified z values
a = np.interp(z, hypso_df['z'], hypso_df['a'])	

def init_mesh(mesh, ic):
	mesh[:, 0] = ic
	
for scenario in scenarios:
	ds = scenarios[scenario]['ds']

	start_date = ds['time'].values[0]
	end_date = ds['time'].values[-1]
	
	sa_bc = ds['so'].values
	theta_bc = ds['thetao'].values
	tis_bc = gsw.t_from_CT(sa_bc, theta_bc, 0)
	
	n_t = ds['time'].size
	
	sa = np.empty((n_z, n_t), order='F')
	tis = np.empty_like(sa)
	
	rho = np.full_like(sa, np.nan)
	kappa = np.full_like(sa, np.nan)
	derived = np.empty((n_z, n_t, bd.indices.n_derived), order='F')

	init_mesh(sa, np.linspace(35.32, 35.24, n_z))
	init_mesh(tis, np.linspace(13.0, 12.2, n_z))
	

	
	# set the time span in seconds
	t_start = 0.0
	t_end = (end_date - start_date).astype('timedelta64[s]').astype(float)
	
	t = np.linspace(t_start, t_end, n_t)
	
	YEAR_SECONDS = 60*60*24*365.25
	
	dt = (t_end - t_start) / (n_t - 1)
	dz = (z[-1] - z[0]) / (n_z - 1)
	
	kappa_buoyancy_coefs = [2e-7, 1.6]

	
	#sa_bc = ds_offshore_doubtful['so'].values
	#tis_bc = gsw.t_from_CT(sa_bc, ds_offshore_doubtful['thetao'].values, 0)
	
	#for scenario in scenarios:
	#	sa_bc = scenarios[scenario]['salinity'].values
	#	tis_bc = gsw.t_from_CT(sa_bc, scenarios[scenario]['temperature'].values, 0)
	#		
	kappa_vals = [5e-4, 1e-3, 8e-3]
	
	output = [{} for _ in range(len(kappa_vals))]

	for i, kappa_val in enumerate(kappa_vals):
		output[i]['sa'] = sa.copy(order='F')
		output[i]['tis'] = tis.copy(order='F')
		output[i]['derived'] = derived.copy(order='F')
	
		bd.basin.compute(z, a, output[i]['sa'], output[i]['tis'], output[i]['derived'], t, sa_bc, tis_bc, bd.computation_modes.normal, bd.diffusivity_modes.constant, bd.renewal_modes.enabled, [kappa_val, np.nan])
	
		sa = output[i]['sa']
		tis = output[i]['tis']
		derived = output[i]['derived']
	
		pres = gsw.p_from_z(z, 0)
	
		sigma0 = gsw.sigma0(sa, tis)
		rho = gsw.rho(sa, tis, np.tile(pres, (n_t, 1)).T)
	
		ds = xr.Dataset(
			{
				'sa': (['depth', 'time'], sa),
				'tis': (['depth', 'time'], tis),
				'rho': (['depth', 'time'], rho),
				'sigma0': (['depth', 'time'], sigma0),
				'pres': (['depth'], pres),
				'renewal': (['depth', 'time'], derived[:, :, bd.indices.renewal-1]),
				'n_freq': (['depth', 'time'], derived[:, :, bd.indices.n_freq-1]),
			},
			coords={
				'depth': np.abs(z),
				'time': ds['time'],
			},
			attrs={
				'kappa': kappa_val,
				'kappa_mode': 'constant'
				}
		)
	
		ds.to_netcdf(f'basin_simulation_{kappa_val}.nc')
	
	## Create a plot of depth averaged properties over time
	fig = plt.figure()
	
	ax_sa = fig.add_subplot(3, 1, 1)
	ax_ct = fig.add_subplot(3, 1, 2)
	ax_sigma0 = fig.add_subplot(3, 1, 3)
	
	sa_mean = np.array([output[i]['sa'].mean(axis=0) for i in range(len(kappa_vals))])
	tis_mean = np.array([output[i]['tis'].mean(axis=0) for i in range(len(kappa_vals))])
	sigma0_mean = np.array([gsw.sigma0(sa_mean[i], tis_mean[i]) for i in range(len(kappa_vals))])
	
	sa_range = [sa_mean.min(axis=0), sa_mean.max(axis=0)]
	tis_range = [tis_mean.min(axis=0), tis_mean.max(axis=0)]
	sigma0_range = [sigma0_mean.min(axis=0), sigma0_mean.max(axis=0)]
	
	ax_sa.plot(ds['time'], sa_mean[1, :], label='$\\kappa ={:.0e}$'.format(kappa_vals[1]), color='blue')
	ax_ct.plot(ds['time'], tis_mean[1, :], label='$\\kappa={:.0e}$'.format(kappa_vals[1]), color='blue')
	ax_sigma0.plot(ds['time'], sigma0_mean[1, :], label='', color='blue')
	
	ax_sa.fill_between(ds['time'], sa_range[0], sa_range[1], color='blue', alpha=0.5)
	ax_ct.fill_between(ds['time'], tis_range[0], tis_range[1], color='blue', alpha=0.5)
	ax_sigma0.fill_between(ds['time'], sigma0_range[0], sigma0_range[1], color='blue', alpha=0.5)
	
	ax_sa.set_ylabel('$S_{ABS}$ \\ $gkg^{-1}$')
	ax_ct.set_ylabel('$T_{CONS}$ \\ $^{\\circ}C$')
	ax_sigma0.set_ylabel('$\\sigma_0$ \\ $kgm^{-3}$')
	
	for ax in [ax_sa, ax_ct, ax_sigma0]:
	    ax.set_xlabel('Time / months')
	    ax.grid('on')
	    ax.legend()
	
	plt.show()
