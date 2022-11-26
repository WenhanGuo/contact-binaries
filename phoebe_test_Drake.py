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

MJD = np.array(df['MJD'])
MJD = MJD - MJD[0]   # set MJD start from 0 for t0 argument
MJD = MJD % 0.3439788   # fold time into delta time

fluxes = 10 ** (-np.array(df['Mag']) + 25)   # obtain flux from mag

# %%
b.set_value('period@binary', value = 0.3439788)
b.set_value('incl', component='binary', value = 89.6)

b['teff@primary'] = 5742
b['teff@secondary'] = 5600

b.flip_constraint('mass@primary', solve_for='sma@binary')
b.set_value('mass@component', component='primary', value=1.25)
b.set_value('q', value = 0.110)

b['requiv@primary'] = 1.37

# %%
print(b.run_checks())   # check if modeling is possible

# %%
b.add_dataset('lc', times=MJD, fluxes=fluxes, dataset='lc01')
b.add_dataset('mesh', compute_times=[0], dataset='mesh01')

b.set_value_all('gravb_bol', 0.32)   # set gravity darkening = 0.32 for both stars since both are convective
b.set_value('pblum_mode', 'dataset-scaled')   # scale passband luminosity to dataset

b.run_compute(model='default')
_ = b.plot(x='phase', show=True)


# %%
b.add_solver('estimator.ebai', ebai_method='knn', solver='ebai_knn')
b.run_solver('ebai_knn', solution='ebai_knn_solution')

# %%
b.flip_constraint('teffratio', solve_for='teff@secondary')
b.flip_constraint('pot@contact_envelope', solve_for='requiv@primary')
# this cell is sus, flipping these might lead model away from true value
print(b.adopt_solution('ebai_knn_solution'))

# %%
b.run_compute(model='ebai_knn_model')
_ = b.plot(x='phase', ls='-', legend=True, show=True)

# %%
