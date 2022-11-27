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
# Reading Drake data for CSS_0113
url = 'https://raw.githubusercontent.com/WenhanGuo/contact-binaries/master/CSS_0113_lcdata.csv'
df = pd.read_csv(url)

MJD = np.array(df['MJD'])
MJD = MJD - MJD[0]   # set MJD start from 0 for t0 argument
MJD = MJD % 0.3439788   # fold time into delta time

fluxes = 10 ** (-np.array(df['Mag']) + 25)   # obtain flux from mag


# %%
b['period@binary'] = 0.3439788   # period = 0.34 day
b['t0_supconj'] = 0.14   # primary eclipse time (zero phase) = 0.14 day
b['incl@binary'] = 89.6

b['teff@primary'] = 5742
b['teff@secondary'] = 5600

b.flip_constraint('mass@primary', solve_for='sma@binary')
b['mass@primary@component'] = 1.25
b['q'] = 0.110

b['requiv@primary'] = 1.37

# %%
print(b.run_checks())   # check if run_compute is possible
print(b)   # check full parameters; it's already close to Christopolou's values

# %%
b.add_dataset('lc', times=MJD, fluxes=fluxes, dataset='lc01')   # add Drake lc dataset
b.add_dataset('orb', compute_times=np.linspace(0,0.3439788,101), dataset='orb01')   # init empty orbit dataset
b.add_dataset('mesh', compute_times=np.linspace(0,0.3439788,31), dataset='mesh01', columns=['teffs'])   # init empty mesh dataset, expose teffs

b.set_value_all('gravb_bol', 0.32)   # set gravity darkening = 0.32 for both stars since both are convective
b.set_value_all('pblum_mode', 'dataset-scaled')   # scale passband luminosity to dataset

b.run_compute(model='default')

# %%
b.plot('lc01', x='phase', size=0.01, legend=True, show=True, save='lc.png')   # plot lc data and forward model
b.plot('mesh01', time=0, legend=True, fc='teffs', ec='None', fcmap='inferno', show=True)   # plot mesh w/ temp color @t0

# animations
b.plot(y={'orb':'ws'}, size=0.01, fc={'mesh':'teffs'}, ec={'mesh':'None'}, 
        fcmap='inferno', animate=True, save='animations_1.gif')   # sync animation for lc, orb, mesh
b.plot('orb01', y='ws', legend=True, animate=True, save='orb2d.gif')   # animate face-on 2d orbit
b.plot('orb01', projection='3d', legend=True, animate=True, save='orb3d.gif')   # animate 3d orbit
b.plot('mesh01', fc='teffs', ec='None', fcmap='inferno', legend=True, animate=True, save='mesh.gif')   # animate mesh

# %%
# start of inverse problem: add and run KNN estimator
b.add_solver('estimator.ebai', ebai_method='knn', solver='ebai_knn')
b.run_solver('ebai_knn', solution='ebai_knn_solution')

# %%
# adopting KNN proposed solutions
# this cell is sus, flipping these might lead model away from true value
b.flip_constraint('teffratio', solve_for='teff@secondary')
b.flip_constraint('pot@contact_envelope', solve_for='requiv@primary')

# print(b.adopt_solution('ebai_knn_solution', trial_run=True))
print(b.adopt_solution('ebai_knn_solution'))

# %%
# forward model from adopted KNN solutions
b.run_compute(model='ebai_knn_model')
_ = b.plot(x='phase', ls='-', legend=True, show=True)

# %%
