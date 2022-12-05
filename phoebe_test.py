# %%
import phoebe
from phoebe import u # units
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# %matplotlib widget
from astropy.timeseries import TimeSeries
from astropy.time import Time

logger = phoebe.logger()

b = phoebe.default_binary(contact_binary=True)

# %%
# Download csv from github, read into pandas
# url = 'https://raw.githubusercontent.com/WenhanGuo/contact-binaries/master/diff_lc.csv'
# df = pd.read_csv(url, delim_whitespace=True)

url = 'https://raw.githubusercontent.com/WenhanGuo/contact-binaries/master/joined_ts_3nights.csv'
df = pd.read_csv(url)
df.set_index(pd.DatetimeIndex(df['time']), inplace=True)
del df['time']

# Convert to astropy TimeSeries
ts = TimeSeries.from_pandas(df)
MJD = ts['time'].mjd
MJD = MJD - MJD[0]   # set MJD start from 0 for t0 argument
MJD = MJD % 0.3439788   # fold time into delta time
fluxes = 10**(-ts['diff_mag']/2.5 + 10)

ts1 = ts.loc['2022-11-06T04:11:20.000':'2022-11-06T09:11:24.000']
MJD1 = ts1['time'].mjd
MJD1 = MJD1 % 0.3439788
fluxes1 = 10**(-ts1['diff_mag']/2.5 + 10)

ts2 = ts.loc['2022-11-21T02:48:28.000':'2022-11-21T08:19:18.000']
MJD2 = ts2['time'].mjd
MJD2 = MJD2 % 0.3439788
fluxes2 = 10**(-ts2['diff_mag']/2.5 + 10)

plt.scatter(MJD1, fluxes1, s=0.5, c='blue')
plt.scatter(MJD2, fluxes2+10**7.72, s=0.5, c='crimson')
# plt.scatter(MJD, fluxes, s=0.5, c='black')
# %%
b.add_dataset('mesh', compute_times=np.linspace(0,0.3439788,31), dataset='mesh01')
b.add_dataset('lc', times=MJD, fluxes=fluxes, dataset='lc01')
b.add_dataset('orb', compute_times=np.linspace(0,0.3439788,101), dataset='orb01')

# %%
print(phoebe.list_online_passbands())
b.set_value('passband', 'SDSS:g')

b.set_value_all('ld_mode', 'lookup')
b.set_value_all('ld_mode_bol', 'lookup')
b.set_value_all('atm', 'ck2004')
b.set_value('pblum_mode', 'dataset-scaled')

b.set_value_all('gravb_bol', 0.32)
b.set_value_all('irrad_frac_refl_bol', 0.5)

b['period@binary'] = 0.3439788   # period = 0.34 day
b['t0_supconj'] = 0.14   # primary eclipse time (zero phase) = 0.14 day
b['incl@binary'] = 89.6
b['Av'] = 0.179

b['teff@primary'] = 5742
b['teff@secondary'] = 5600

b.flip_constraint('mass@primary', solve_for='sma@binary')
b['mass@primary@component'] = 1.25
b['q'] = 0.110

b['requiv@primary'] = 1.37

print(b)

# %%
# b.run_compute(model='default')
_ = b.plot(x='times', show=True, size=0.0015)

# %%
b.add_solver('estimator.lc_periodogram')
b.run_solver(kind='lc_periodogram', lc_datasets='lc01')

# %%
b.add_solver('estimator.ebai', ebai_method='knn', solver='ebai_knn', overwrite=True)
b.run_solver('ebai_knn', solution='ebai_knn_solution', phase_bin=False)

# %%
b.flip_constraint('teffratio', solve_for='teff@secondary')
b.flip_constraint('pot@contact_envelope', solve_for='requiv@primary')

print(b.adopt_solution('ebai_knn_solution'))
# print(b.adopt_solution('ebai_knn_sol', adopt_parameters=['t0_supconj','teffratio','incl']))

# %%
b.run_compute(model='ebai_knn_model', overwrite=True)
_ = b.plot('lc01', x='phase', ls='-', legend=True, show=True)

# %%
b.add_solver('optimizer.nelder_mead', 
            fit_parameters=['teffratio', 'incl@binary', 'q', 'per0'], solver='nm_solver')
b.run_solver('nm_solver', maxiter=10000, solution='nm_sol')

# %%
b.add_solver('sampler.emcee', solver='emcee_solver')
b.set_value('compute', solver='emcee_solver', value='fastcompute')
b.set_value('pblum_mode', 'dataset-coupled')

b.add_distribution({'t0_supconj': phoebe.gaussian_around(0.01),
                    'teffratio@binary': phoebe.gaussian_around(0.1),
                    'incl@binary': phoebe.gaussian_around(5),
                    'fillout_factor@contact_envelope': phoebe.gaussian_around(0.5),
                    'q@primary': phoebe.gaussian_around(0.5),
                    'pblum@primary': phoebe.gaussian_around(0.2),
                    'sigmas_lnf@lc01': phoebe.uniform(-1e9, -1e4),
                   }, distribution='ball_around_guess')
b.run_compute(compute='fastcompute', sample_from='ball_around_guess',
                sample_num=10, model='init_from_model')
_ = b.plot('lc01', x='phase', ls='-', model='init_from_model', show=True)

# %%
b['init_from'] = 'ball_around_guess'
b.set_value('nwalkers', solver='emcee_solver', value=14)
b.set_value('niters', solver='emcee_solver', value=500)

b.run_solver('emcee_solver', solution='emcee_solution')
print(b.adopt_solution(solution='emcee_solution', distribution='emcee_posteriors'))
_ = b.plot_distribution_collection(distribution='emcee_posteriors', show=True)
b.uncertainties_from_distribution_collection(distribution='emcee_posteriors', sigma=3, tex=False)

# %%
_ = b.plot(solution='emcee_solution', style='lnprobability',
            burnin=100, thin=1, lnprob_cutoff=3600,
            show=True)

# We can fix the following if we know which values result in a nice chain
b.set_value('burnin', 100)
b.set_value('thin', 1)
b.set_value('lnprob_cutoff', 3600)

# %%
# Show corner plot
_ = b.plot(solution='emcee_solution', style='corner', show=True)

# %%
# Show corner plot with failed and rejected samples
_ = b.plot(solution='emcee_solution', style='failed', show=True)

# %%
# Show history of each sampled parameter for all walkers
_ = b.plot(solution='emcee_solution', style='trace', show=True)

# %%
b.run_compute(compute='fastcompute', sample_from='emcee_solution',
                sample_num=20, model='emcee_sol_model')
_ = b.plot('lc01', x='phase', ls='-', model='emcee_sol_model', show=True)