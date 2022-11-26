##  http://phoebe-project.org/docs/2.4/physics

## Glossary
vgamma          systemic velocity (similar to radial velocity)
ltte            correct for light travel time effect [True/False]
t0              time at which all values are provided and all computations start
l3_mode         third light mode: flux unit or fraction of total flux ['flux'/'fraction']
    l3              appear if l3_mode == 'flux'. Third light in flux units
    l3_frac         appear if l3_mode == 'fraction'. Third light as a fraction of total flux
distance        distance between system center and observer at t0
ebv             color excess E(B - V)
Av              extinction Av
Rv              extinction law parameter (default = 3.1)

t0_supconj      time of primary eclipse (superior conjunction), usually defined as the zero-phase
t0_perpass      time of closest approach between components (periastron passage)
t0_ref          time at which the primary passes an arbitrary reference point (same as HJD0)
ecc             eccentricity
ecosw           eccentricity times cos of periastron angle
esinw           eccentricity times sin of periastron angle
period          sidereal orbital period, zero point defined near one of the t0s
dperdt          dperdt = 0 if not precessing, dperdt ≠ 0 if precessing (apsidal motion)
                dperdt between t0 and other reference times is used in constrainining t0_subconj, t0_perpass, t0_ref
period_anom     anomalistic orbital period (time between two successive periastron passages)
                longer than sidereal period if dperdt ≠ 0 --> breaks Kepler's law assumptions
per0            periastron angle (angle between periastron and ascending node), defined at t0
syncpar         ratio between sidereal orbital period and rotational period wrt the sky
                syncpar = 1 and dperdt = 0 --> stars are co-rotating
                syncpar = 1 and dperdt ≠ 0 --> star is rotating wrt the sky at the same rate as sidereal period
dpdt            time derivative of anomalistic period (= sidereal period if not precessing)
                set dpdt = 0 to see eclipses spread across the phase-space when plotting
phases_dpdt     dpdt used in conversion between compute_times and compute_phases
                when mapping, set phases_dpdt = 'none' so that compute_times are direct multiples of period
pitch           misalignment of a star in the direction of inclination
yaw             misalignment in the direction of longitude of a star's equator

requiv          equivalent radius
pot             potential of the contact envelope
fillout_factor  fillout_factor of the envelope
logg            stellar surface gravity at requiv
spot parameters
    colat           latitude of spot, 0 is defined as the North (spin) Pole
    long            longitude of spot, 0 is defined as pointing towards the other star
    radius          angular radius of spot
    relteff         ratio of temperature of the spot to the local intrinsic value
    enabled         if spot is enabled [True/False]
eclipse_method  method of eclipse calculation ['native'/'visible_partial']
                native: computes what percentage (by area) of each triangle is visible
                visible_partial: assigns visibility = 0.5 to partially visible triangles
pblum_mode      mode to handle passband luminosity
    component-coupled   provide pblum for one star (by default L1), compute pblum for the other
    decoupled           provide pblums for each star independently
    absolute            obtain unscaled pblums, in passband watts, computed from atmosphere tables
    dataset-scaled      calculate pblum for each star from abs flux, scale to dataset
    dataset-coupled     all datasets scaled with the same scaling factor
gravb_bol       the β coefficient for gravity darkening corrections
irrad_frac_refl_bol     fraction of incident flux handled by reflection/irradiation (heating w/o redistribution)
irrad_frac_lost_bol     fraction of incident flux lost/ignored
irrad_method    the method to use to handle all irradiation effects ['none'/'wilson'/'horvat']
rv_offset       radial velocity offset

atm             per-component atmosphere table ['ck2004'/'blackbody'/'extern_atmx'/'extern_planckint'/'phoenix']
                if not 'ck2004' (Castelli-Kurucz): ld_func must not be 'interp'
passband        passband to compute intensities for each star
ld_mode         mode to use for per-component passband limb darkening ['interp'/'lookup'/'manual']
ld_mode_bol     mode to use for per-component bolometric limb darkening ['lookup'/'manual']
ld_func_bol     per-component bolometric limb darkening model
                ['linear'/'logarithmic'/'quadratic'/'square_root'/'power']
ld_coeffs_source_bol    source for bolometric limb darkening coefficients ['auto'/'phoenix'/'ck2004']
ld_coeffs_bol   per-component bolometric limb darkening coefficients
intens_weighting        whether passband intensities are weighted by energy or photons ['energy'/'photon']
exptime         exposure time (time is defined as mid-exposure)
fti_method      how to handle finite-time integration (when non-zero exptime) ['none'/'oversample']
fti_oversample  number of times to sample per-datapoint (averaged) for finite-time integration, default to 5
## ------------------------------------- 1. System Effects --------------------------------------
Systemic Velocity [vgamma] [ltte] [t0]
    <!-- http://phoebe-project.org/docs/2.4/tutorials/vgamma -->
Rømer & Light Travel Time Effects [ltte]
    <!-- http://phoebe-project.org/docs/2.4/tutorials/ltte -->
Third Light [l3_mode] [l3] [l3_frac]
    <!-- http://phoebe-project.org/docs/2.4/tutorials/l3 -->
    ** 3rd light is only a fixed flux offset (extraneous light added to the system). NOT a 3rd body
Distance [distance]
    <!-- http://phoebe-project.org/docs/2.4/tutorials/distance -->
Extinction [ebv] [Av] [Rv]
    <!-- http://phoebe-project.org/docs/2.4/tutorials/ebv_Av_Rv -->
    2F, 1C
    default C: [ebv]

## ------------------------------------- 2. Orbital Effects -------------------------------------
Various t0s [t0_supconj] [t0_perpass] [t0_ref] [t0]
    <!-- http://phoebe-project.org/docs/2.4/tutorials/t0s -->
    1F, 2C
    default F: [t0_supconj]
Eccentricity & Volume Conservation [ecc] [ecosw] [esinw]
    <!-- http://phoebe-project.org/docs/2.4/tutorials/ecc -->
    1F, 2C
    default F: [ecc]
Apsidal Motion [period] [dperdt] [period_anom] [per0] [syncpar]
    <!-- http://phoebe-project.org/docs/2.4/tutorials/apsidal_motion -->
Period Change [dpdt]
    <!-- http://phoebe-project.org/docs/2.4/tutorials/dpdt -->
Misalignment [pitch] [yaw]
    <!-- http://phoebe-project.org/docs/2.4/tutorials/pitch_yaw -->
    2F, 1C
    default F: orbital inclination, pitch / orbital long_an, yaw
    default C: inclination of component / long_an of component

## ------------------------------------- 3. Stellar Effects -------------------------------------
Equivalent Radius [requiv]
    <!-- http://phoebe-project.org/docs/2.4/tutorials/requiv -->
Critical Radii: Contact Systems [requiv@primary] [requiv@secondary] [pot] [fillout_factor]
    <!-- http://phoebe-project.org/docs/2.4/tutorials/requiv_crit_contact -->
    1F, 3C
    default F: [requiv@primary]
    default C: [requiv@secondary] [pot] [fillout_factor]
    [requiv] and [pot] are also constrained by C params [requiv_max] [requiv_min] [pot_max] [pot_min]
    To provide any of the three C params, see link to flip constraints
Surface Gravity [logg]
    <!-- http://phoebe-project.org/docs/2.4/tutorials/logg -->
Spots [colat] [long] [radius] [relteff] [enabled]
    <!-- http://phoebe-project.org/docs/2.4/tutorials/spots -->
Eclipse Detection [eclipse_method]
    <!-- http://phoebe-project.org/docs/2.4/tutorials/eclipse -->
Passband Luminosity [pblum_mode] [pblum]
    <!-- http://phoebe-project.org/docs/2.4/tutorials/pblum -->
    [pblum_mode] = 'component-coupled' (default)
        User provides passband luminosities for a single star in the system, secondary star is scaled accordingly.
        default: [pblum] = 4π (12.57 W). (flux from primary of ~1; identical secondary so total flux ~2)
        ** setting component temperature has a large effect since flux scales to T^4
    [pblum_mode] = 'decoupled' (undesirable)
        Separate [pblum] for primary and secondary; both params are available now and can have different values.
        ** setting component temperature has no effect since luminosity will be rescaled to primary
    [pblum_mode] = 'absolute'
        Luminosities and fluxes will be returned in absolute units and not rescaled.
    [pblum_mode] = 'dataset-scaled' (may be desired)
        Only allowed if fluxes are attached to the dataset. Resulting model will be scaled to best fit the data.
    [pblum_mode] = 'dataset-coupled' (may be desired)
        Allows for the same scaling factor to be applied to two different datasets. 
Gravity Brightening/Darkening [gravb_bol]
    <!-- http://phoebe-project.org/docs/2.4/tutorials/gravb_bol -->
    teff 8000+ (radiative atm): [gravb_bol] >= 0.9 (suggest 1.0)
    teff 6600-8000 (intermittent): [gravb_bol] 0.32-1.0
    teff 6600- (convective atm): [grav_bol] < 0.9 (suggest 0.32)
Reflection & Heating [irrad_frac_refl_bol] [irrad_frac_lost_bol]
    <!-- http://phoebe-project.org/docs/2.4/tutorials/reflection_heating -->
    Dictates what fraction of incident light is handled by reflection / lost flux
    for each component, [irrad_frac_refl_bol] + [irrad_frac_lost_bol] = 1
    1F, 1C
    default F: [irrad_frac_refl_bol]
Reflection & Heating: Lambert Scattering [irrad_method]
    <!-- http://phoebe-project.org/docs/2.4/tutorials/irrad_method_horvat -->
    [irrad_method] = 'none': ignore irradiation
    [irrad_method] = 'wilson': Wilson's original reflection scheme
    [irrad_method] = 'horvat': the new Horvat scheme which includes Lambert Scattering
Radial Velocity Offsets [rv_offset]
    <!-- http://phoebe-project.org/docs/2.4/tutorials/rv_offset -->

## --------------------------- 4. Passband/Atmosphere/Dataset Effects ---------------------------
Passbands & Atmospheres [atm] [passband]
    <!-- http://phoebe-project.org/docs/2.4/tutorials/atm_passbands -->
Limb Darkening [ld_mode] [ld_mode_bol] [ld_func] [ld_func_bol] [ld_coeffs_source] [ld_coeffs_source_bol] [ld_coeffs] [ld_coeffs_bol]
    <!-- http://phoebe-project.org/docs/2.4/tutorials/limb_darkening -->
    [ld_mode_bol] = 'lookup': interpolate bolometric LD coefficients per-component
    [ld_mode_bol] = 'manual': pass bolometric LD coefficients manually
    [ld_coeffs_source_bol] = 'auto': interpolate from applicable table according to atm
    [ld_mode] = 'interp'
        Interpolate passband LD coefficients from atmospheric tables
        [ld_func], [ld_coeffs_source] and [ld_coeffs] are invisible (irrelevant)
    [ld_mode] = 'lookup': interpolate passband LD coefficients per-element, expose [ld_func] and [ld_coeffs_source]
    [ld_mode] = 'manual': manually set [ld_coeffs], expose [ld_func] and [ld_coeffs], need to fit LD model
Gravitational Redshift (irrelevant for lc) [rv_grav]
    <!-- http://phoebe-project.org/docs/2.4/tutorials/grav_redshift -->
Intensity Weighting [intens_weighting]
    <!-- http://phoebe-project.org/docs/2.4/tutorials/intens_weighting -->
Finite Time of Integration (fti) [exptime] [fti_method]
    <!-- http://phoebe-project.org/docs/2.4/tutorials/fti -->
    [fti_method] = 'oversample': expose [fti_oversample]
Gaussian Processes
    <!-- http://phoebe-project.org/docs/2.4/examples/minimal_GPs -->