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
b = b


# %%
directory = '/Users/danny/My_Root/ASTRO/Contact_Binary/data/CSS_034852/20220719/diff_lc.ecsv'
ts = TimeSeries.read(directory, time_column='time')

# %%
b.add_dataset('mesh', compute_times=[0], dataset='mesh01')
b.add_dataset('lc', times=ts.time.jd, fluxes=ts['diff_flux'], dataset='lc01')

b.set_value_all('ld_mode', 'manual')
b.set_value_all('ld_mode_bol', 'manual')
b.set_value_all('atm', 'blackbody')
b.set_value('pblum_mode', 'dataset-scaled')

b.run_compute(model='default')
_ = b.plot(x='phase', show=True)

# %%
b.add_solver('estimator.ebai', ebai_method='knn', solver='ebai_knn')
b.run_solver('ebai_knn', solution='ebai_knn_solution')

# %%
# b.flip_constraint('teffratio', solve_for='teff@secondary')
# b.flip_constraint('pot@contact_envelope', solve_for='requiv@primary')

print(b.adopt_solution('ebai_knn_solution'))

# %%
b.run_compute(model='ebai_knn_model')
_ = b.plot(x='phase', ls='-', legend=True, show=True)

# %%
