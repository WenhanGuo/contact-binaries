# %%
import phoebe
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

logger = phoebe.logger()
b = phoebe.default_binary(contact_binary=True)

# %%
# Reading Drake data for CSS_0113
url = 'https://raw.githubusercontent.com/WenhanGuo/contact-binaries/master/CSS_0113_lcdata.csv'
df = pd.read_csv(url)

MJD = np.array(df['MJD'])
MJD = MJD - MJD[0]   # set MJD start from 0 for t0 argument
MJD = MJD % 0.3439788   # fold time into delta time

fluxes = 10 ** (-np.array(df['Mag'])/2.5 + 10)   # obtain flux from mag

# %%
orbphases = phoebe.linspace(0,1,101)
meshphases = phoebe.linspace(0,1,31)
b.add_dataset('lc', times=MJD, fluxes=fluxes, dataset='lc01')   # add Drake lc dataset
b.add_dataset('orb', compute_phases=orbphases, dataset='orb01')   # init empty orbit dataset
b.add_dataset('mesh', compute_phases=meshphases, dataset='mesh01', columns=['teffs'])   # init empty mesh dataset, expose teffs

# %%
b.set_value_all('pblum_mode', 'dataset-scaled')   # scale passband luminosity to dataset
b.set_value_all('gravb_bol', 0.32)   # set gravity darkening = 0.32 for both stars since both are convective
b.set_value_all('irrad_frac_refl_bol', 0.5)   # set surface albedo = 0.5

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

print(b.run_checks())   # check if run_compute is possible
print(b)   # check full parameters; it's already close to Christopolou's values

# %%
b.run_compute(model='default')

# %%
# simple plotting
b.plot('lc01', x='phase', size=0.012, legend=True, show=True, save='./cb_visu_drake/lc.png')   # plot lc data and forward model
b.plot('mesh01', phase=0, legend=True, fc='teffs', ec='None', fcmap='viridis', show=True)   # plot mesh w/ temp color @t0
# animations
b.plot(y={'orb':'ws'}, size=0.01, fc={'mesh':'teffs'}, ec={'mesh':'None'}, 
        fcmap='viridis', animate=True, save='./cb_visu_drake/animations_sync.gif')   # sync animation for lc, orb, mesh
b.plot('orb01', y='ws', legend=True, animate=True, save='./cb_visu_drake/orb2d.gif')   # animate face-on 2d orbit
b.plot('orb01', projection='3d', legend=True, animate=True, save='./cb_visu_drake/orb3d.gif')   # animate 3d orbit
b.plot('mesh01', fc='teffs', ec='None', fcmap='viridis', legend=True, animate=True, save='./cb_visu_drake/mesh.gif')   # animate mesh

# %%
# b.add_solver('estimator.lc_periodogram')
# b.run_solver(kind='lc_periodogram', lc_datasets='lc01')

# %%
# start of inverse problem: add and run KNN estimator
b.add_solver('estimator.ebai', ebai_method='knn', solver='ebai_knn')
b.run_solver('ebai_knn', solution='ebai_knn_sol')
print(b.adopt_solution('ebai_knn_sol', trial_run=True))   # see proposed KNN solution params before adopting

# %%
b.flip_constraint('teffratio', solve_for='teff@secondary')

# if adopt all proposed params, uncomment below:
b.flip_constraint('pot@contact_envelope', solve_for='requiv@primary')
print(b.adopt_solution('ebai_knn_solution'))

# if not adopting q, uncomment below:
# print(b.adopt_solution('ebai_knn_sol', adopt_parameters=['t0_supconj','teffratio','incl']))

# %%
# forward model from adopted KNN solutions
b.add_dataset('mesh', compute_phases=meshphases, dataset='mesh02', columns=['teffs'])
b.run_compute(model='ebai_knn_model', overwrite=True)
b.plot('lc01', x='phase', ls='-', legend=True, show=True)
b.plot('mesh02', fc='teffs', ec='None', fcmap='viridis', legend=True, animate=True, save='./cb_visu_drake/mesh_inverse_drake.gif')

# %%
# nelder_mead optimizer run with 1000 iterations
b.add_solver('optimizer.nelder_mead', 
            fit_parameters=['teffratio', 'incl@binary', 'q', 'per0'], solver='nm_solver')
b.run_solver('nm_solver', maxiter=10000, solution='nm_sol')
