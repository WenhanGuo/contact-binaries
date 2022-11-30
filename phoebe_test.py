# %%
import phoebe
from phoebe import u # units
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# %matplotlib widget
from astropy.timeseries import TimeSeries

logger = phoebe.logger()

b = phoebe.default_binary(contact_binary=True)


# %%
# Download csv from github, read into pandas
# url = 'https://raw.githubusercontent.com/WenhanGuo/contact-binaries/master/diff_lc.csv'
# df = pd.read_csv(url, delim_whitespace=True)

url = 'https://raw.githubusercontent.com/WenhanGuo/contact-binaries/master/joined_ts_2nights.csv'
df = pd.read_csv(url)
df.set_index(pd.DatetimeIndex(df['time']), inplace=True)
del df['time']

# Convert to astropy TimeSeries
ts = TimeSeries.from_pandas(df)

MJD = ts['time'].mjd
MJD = MJD - MJD[0]   # set MJD start from 0 for t0 argument
MJD = MJD % 0.3439788   # fold time into delta time

fluxes = 10**(-ts['diff_mag']/2.5 + 10)

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
#b.run_compute(model='default')
_ = b.plot(x='times', show=True, size=0.001)

# %%
b.add_solver('estimator.lc_periodogram')
b.run_solver(kind='lc_periodogram', lc_datasets='lc01')

# %%
b.add_solver('estimator.ebai', ebai_method='knn', solver='ebai_knn', overwrite=True)
b.run_solver('ebai_knn', solution='ebai_knn_solution', phase_bin=False)

# %%
# b.flip_constraint('teffratio', solve_for='teff@secondary')
# b.flip_constraint('pot@contact_envelope', solve_for='requiv@primary')

print(b.adopt_solution('ebai_knn_solution'))
# print(b.adopt_solution('ebai_knn_sol', adopt_parameters=['t0_supconj','teffratio','incl']))

# %%
b.run_compute(model='ebai_knn_model')
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
                    'fillout_factor@contact_envelope': phoebe.gaussian_around(0.4),
                    'q@primary': phoebe.gaussian_around(0.2),
                    'pblum@primary': phoebe.gaussian_around(0.2),
                    'sigmas_lnf@lc01': phoebe.uniform(-1e9, -1e4),
                   }, distribution='ball_around_guess')