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
# setting params
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
# forward model
orbphases = phoebe.linspace(0,1,101)
meshphases = phoebe.linspace(0,1,31)
b.add_dataset('lc', times=MJD, fluxes=fluxes, dataset='lc01')   # add Drake lc dataset
b.add_dataset('orb', compute_phases=orbphases, dataset='orb01')   # init empty orbit dataset
b.add_dataset('mesh', compute_phases=meshphases, dataset='mesh01', columns=['teffs'])   # init empty mesh dataset, expose teffs

b.set_value_all('gravb_bol', 0.32)   # set gravity darkening = 0.32 for both stars since both are convective
b.set_value_all('irrad_frac_refl_bol', 0.5)   # set surface albedo = 0.5
b.set_value_all('pblum_mode', 'dataset-scaled')   # scale passband luminosity to dataset

b.run_compute(model='default')

# %%
# simple plotting
b.plot('lc01', x='phase', size=0.01, legend=True, show=True, save='lc.png')   # plot lc data and forward model
b.plot('mesh01', phase=0, legend=True, fc='teffs', ec='None', fcmap='inferno', show=True)   # plot mesh w/ temp color @t0
# animations
b.plot(y={'orb':'ws'}, size=0.01, fc={'mesh':'teffs'}, ec={'mesh':'None'}, 
        fcmap='inferno', animate=True, save='animations_1.gif')   # sync animation for lc, orb, mesh
b.plot('orb01', y='ws', legend=True, animate=True, save='orb2d.gif')   # animate face-on 2d orbit
b.plot('orb01', projection='3d', legend=True, animate=True, save='orb3d.gif')   # animate 3d orbit
b.plot('mesh01', fc='teffs', ec='None', fcmap='inferno', legend=True, animate=True, save='mesh.gif')   # animate mesh

# %%
# start of inverse problem: add and run KNN estimator
b.add_solver('estimator.ebai', ebai_method='knn', solver='ebai_knn')
b.run_solver('ebai_knn', solution='ebai_knn_sol')
print(b.adopt_solution('ebai_knn_sol', trial_run=True))   # see proposed KNN solution params before adopting

# %%
# adopting KNN proposed solutions
b.flip_constraint('teffratio', solve_for='teff@secondary')
print(b.adopt_solution('ebai_knn_sol', adopt_parameters=['t0_supconj','teffratio','incl']))

# flipping pot and adopting all params might lead model away from true value
# b.flip_constraint('pot@contact_envelope', solve_for='requiv@primary')
# print(b.adopt_solution('ebai_knn_sol'))

# %%
# forward model from adopted KNN solutions
b.run_compute(model='ebai_knn_model')
b.plot('lc01', x='phase', ls='-', legend=True, show=True)

# %%
# nelder_mead optimizer run with 1000 iterations
b.add_solver('optimizer.nelder_mead', 
            fit_parameters=['teffratio', 'incl@binary', 'q', 'per0'], solver='nm_solver')
b.run_solver('nm_solver', maxiter=1000, solution='nm_sol')

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
# Show corner plot
_ = b.plot(solution='emcee_solution', style='corner', show=True)

# %%
# Show corner plot with failed and rejected samples
_ = b.plot(solution='emcee_solution', style='failed', show=True)

# %%
b.run_compute(compute='fastcompute', sample_from='emcee_solution',
                sample_num=20, model='emcee_sol_model')
_ = b.plot('lc01', x='phase', ls='-', model='emcee_sol_model', show=True)