##  http://phoebe-project.org/docs/2.4/physics

## Glossary
vgamma      systemic velocity (similar to radial velocity)
ltte        correct for light travel time effect [True/False]
t0          time at which all values are provided and all computations start
l3_mode     third light mode: flux unit or fraction of total flux ['flux'/'fraction']
    l3          appear if l3_mode == 'flux'. Third light in flux units
    l3_frac     appear if l3_mode == 'fraction'. Third light as a fraction of total flux
distance    distance between system center and observer at t0
ebv         color excess E(B - V)
Av          extinction Av
Rv          extinction law parameter (default = 3.1)

t0_supconj  time of primary eclipse (superior conjunction), usually defined as the zero-phase
t0_perpass  time of closest approach between components (periastron passage)
t0_ref      time at which the primary passes an arbitrary reference point (same as HJD0)
ecc         eccentricity
ecosw       eccentricity times cos of periastron angle
esinw       eccentricity times sin of periastron angle
period      sidereal orbital period, zero point defined near one of the t0s
dperdt      dperdt = 0 if not precessing, dperdt ≠ 0 if precessing (apsidal motion)
            dperdt between t0 and other reference times is used in constrainining t0_subconj, t0_perpass, t0_ref
period_anom anomalistic orbital period (time between two successive periastron passages)
            longer than sidereal period if dperdt ≠ 0 --> breaks Kepler's law assumptions
per0        periastron angle (angle between periastron and ascending node), defined at t0
syncpar     ratio between sidereal orbital period and rotational period wrt the sky
            syncpar = 1 and dperdt = 0 --> stars are co-rotating
            syncpar = 1 and dperdt ≠ 0 --> star is rotating wrt the sky at the same rate as sidereal period
dpdt        time derivative of anomalistic period (= sidereal period if not precessing)
            set dpdt = 0 to see eclipses spread across the phase-space when plotting
phases_dpdt dpdt used in conversion between compute_times and compute_phases
            when mapping, set phases_dpdt = 'none' so that compute_times are direct multiples of period
pitch       misalignment of a star in the direction of inclination
yaw         misalignment in the direction of longitude of a star's equator

requiv      equivalent radius
pot         potential of the contact envelope
fillout_factor  fillout_factor of the envelope
logg        stellar surface gravity at requiv
spot params
    colat   latitude of spot, 0 is defined as the North (spin) Pole
    long    longitude of spot, 0 is defined as pointing towards the other star
    radius  angular radius of spot
    relteff ratio of temperature of the spot to the local intrinsic value
    enabled if spot is enabled [True/False]
eclipse_method  method of eclipse calculation ['native'/'visible_partial']
                native: computes what percentage (by area) of each triangle is visible
                visible_partial: assigns visibility = 0.5 to partially visible triangles


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
    <!-- http://phoebe-project.org/docs/2.4/tutorials/logghttp://phoebe-project.org/docs/2.4/tutorials/logg -->
Spots [colat] [long] [radius] [relteff] [enabled]
    <!-- http://phoebe-project.org/docs/2.4/tutorials/spots -->
Eclipse Detection [eclipse_method]
    <!-- http://phoebe-project.org/docs/2.4/tutorials/eclipse -->
Passband Luminosity
    <!-- http://phoebe-project.org/docs/2.4/tutorials/pblum -->
Gravity Brightening/Darkening
Reflection & Heating
Reflection & Heating: Lambert Scattering
Radial Velocity Offsets

## --------------------------- 4. Passband/Atmosphere/Dataset Effects ---------------------------
Passbands & Atmospheres
Passband Luminosity
Limb Darkening
Third Light
Gravitational Redshift
Radial Velocity Offsets
Intensity Weighting
Finite Time of Integration
Gaussian Processes