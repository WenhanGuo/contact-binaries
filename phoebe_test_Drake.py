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
# This cell is for reading Drake data for CSS_0113
url = 'https://raw.githubusercontent.com/WenhanGuo/contact-binaries/master/CSS_0113_lcdata.csv'
df = pd.read_csv(url)


# %%
b.set_value('period@binary', value = 0.3439788)
b.set_value('teff@component', component='primary', value=5742.)
b.set_value('teff@component', component='secondary', value=5600.)
b.set_value('incl', component='binary', value = 89.6)
# b.set_value('q', value = 0.110)
# b.set_value('requiv@component', component='primary', value=0.594)
# b.set_value('requiv@component', component='secondary', value=0.237)

# b.set_value('fillout_factor@component@envelope', value=0.68)
# b.flip_constraint('mass@primary', solve_for='sma')
# b.set_value('mass@component', component='primary', value=1.25)
# b.set_value('mass@component', component='secondary', value=0.14)

# %%
MJD = np.array(df['MJD'])
fluxes = 10 ** (-np.array(df['Mag']) + 25)
fluxes = fluxes / np.median(fluxes)

b.add_dataset('lc', times=MJD, fluxes=fluxes, dataset='lc01')
b.add_dataset('mesh', compute_times=[0], dataset='mesh01')

b.set_value_all('ld_mode', 'manual')
b.set_value_all('ld_mode_bol', 'manual')
b.set_value_all('atm', 'blackbody')
# b.set_value('pblum_mode', 'dataset-scaled')

b.run_compute(model='default')
_ = b.plot(x='phase', show=True)


# %%
b.add_solver('estimator.ebai', ebai_method='knn', solver='ebai_knn')
b.run_solver('ebai_knn', solution='ebai_knn_solution')

# %%
b.flip_constraint('teffratio', solve_for='teff@secondary')
b.flip_constraint('pot@contact_envelope', solve_for='requiv@primary')

print(b.adopt_solution('ebai_knn_solution'))

# %%
b.run_compute(model='ebai_knn_model')
_ = b.plot(x='phase', ls='-', legend=True, show=True)

# %%
