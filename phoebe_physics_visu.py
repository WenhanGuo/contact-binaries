# %%
import phoebe
import matplotlib.pyplot as plt
import numpy as np

logger = phoebe.logger()
b = phoebe.default_star()

# %%
# lcphases and meshphases settings (global). Changing these changes frame number in animation
lcphases = phoebe.linspace(0,1,51)
meshphases = phoebe.linspace(0,1,51)
b.add_dataset('lc', compute_phases=lcphases, dataset='lc01')
b.add_dataset('mesh', compute_phases=meshphases, dataset='mesh01', columns=['teffs'])

# %%
# single star bare mesh, no color
b.set_value_all('gravb_bol', 0)
b.run_compute(model='default', overwrite=True)
b.plot('mesh01', phase=0, show=True, save='./physics_visu/1star_bare.png')

# %%
# single star mesh, gravity darkening & teffs color
b.set_value_all('gravb_bol', 0.32)
b.run_compute(model='default', overwrite=True)
b.plot('mesh01', phase=0, fc='teffs', ec='None', fcmap='viridis', show=True, save='./physics_visu/1star_gd.png')

# %%
# single star bare mesh with spot, animation
b.add_spot(radius=30, colat=80, long=0, relteff=0.9)
b.run_compute(model='default', distortion_method='rotstar', overwrite=True)
b.plot('mesh01', fc='teffs', animate=True, save='./physics_visu/1star_spot.gif')

# %%
# binary bare mesh animation: 1. no roche, 2. with roche 
b = phoebe.default_binary()
b['requiv@primary@component'] = 2
b.add_dataset('mesh', compute_phases=meshphases, dataset='mesh01', columns=['teffs'])
b.run_compute(distortion_method='roche', model='rochemodel')
b.run_compute(distortion_method='rotstar', model='rotstarmodel')
b.plot('mesh01', model='rotstarmodel', animate=True, save='./physics_visu/binary_noroche.gif')
b.plot('mesh01', model='rochemodel', animate=True, save='./physics_visu/binary_roche.gif')

# %%
# binary (almost contact) animation: 1. bare mesh, 2. teffs color, 3. with lc
b = phoebe.default_binary()
b['requiv@primary@component'] = 1.9
b['requiv@secondary@component'] = 2
b['teff@primary'] = 6000
b['teff@secondary'] = 5800
b.add_dataset('lc', compute_phases=lcphases, dataset='lc01')
b.add_dataset('mesh', compute_phases=meshphases, dataset='mesh01', columns=['teffs'])
b.run_compute(model='default', overwrite=True)
b.plot('mesh01', model='default', animate=True, save='./physics_visu/binary_close.gif')
b.plot('mesh01', model='default', fc='teffs', ec='None', animate=True, fcmap='viridis', save='./physics_visu/binary_close_temp.gif')
b.plot(model='default', ylim={'lc':(0.75,2.1)}, phases=meshphases, size=0.012, fc={'mesh':'teffs'}, ec={'mesh':'None'}, 
        fcmap='viridis', animate=True, save='./physics_visu/binary_close_lc.gif')

# %%
# binary (almost contact) animation: 4. with inclination
b['orbit@component@incl'] = 70
b.run_compute(model='incl70', overwrite=True)
b.plot(model='incl70', ylim={'lc':(0.75,2.1)}, phases=meshphases, size=0.012, fc={'mesh':'teffs'}, ec={'mesh':'None'}, 
        fcmap='viridis', animate=True, save='./physics_visu/binary_incl70.gif')

# %%
# contact binary animation: equal masses
b = phoebe.default_binary(contact_binary=True)
b.add_dataset('lc', compute_phases=lcphases, dataset='lc01')
b.add_dataset('mesh', compute_phases=meshphases, dataset='mesh01', columns=['teffs'])
b.run_compute(model='default', overwrite=True)
b.plot(phases=meshphases, size=0.012, fc={'mesh':'teffs'}, ec={'mesh':'None'}, 
        fcmap='viridis', animate=True, save='./physics_visu/contact_default.gif')

# %%
# contact binary animation: low mass ratio
b = phoebe.default_binary(contact_binary=True)
b['period@binary'] = 0.3439788
b.flip_constraint('mass@primary', solve_for='sma@binary')
b['mass@primary@component'] = 1.25
b['q'] = 0.110
b['requiv@primary'] = 1.37
b.add_dataset('lc', compute_phases=lcphases, dataset='lc01')
b.add_dataset('mesh', compute_phases=meshphases, dataset='mesh01', columns=['teffs'])
b.run_compute(model='default', overwrite=True)
b.plot(phases=meshphases, size=0.012, fc={'mesh':'teffs'}, ec={'mesh':'None'}, 
        fcmap='viridis', animate=True, save='./physics_visu/contact_lowq.gif')

# %%
