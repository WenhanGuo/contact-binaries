##  http://phoebe-project.org/docs/2.4/physics

## Glossary
vgamma      systemic velocity (similar to radial velocity)
ltte        correct for light travel time effects [True/False]
t0          time at which all values are provided, should set to the start of time data
l3_mode     third light mode: flux unit or fraction of total flux ['flux'/'fraction']
    l3          appear if l3_mode == 'flux'. Third light in flux units
    l3_frac     appear if l3_mode == 'fraction'. Third light as a fraction of total flux
distance    distance between system center and observer at t0
ebv         extinction E(B - V)
Av          extinction Av
Rv          extinction law parameter (default = 3.1)


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
Various t0s
Eccentricity & Volume Conservation
Apsidal Motion
Period Change
Misalignment
Rømer & Light Travel Time Effects
Beaming & Boosting

## ------------------------------------- 3. Stellar Effects -------------------------------------
Equivalent Radius
Potentials
Critical Radii: Detached Systems
Critical Radii: Semidetached Systems
Critical Radii: Contact Systems [requiv@primary] [requiv@secondary] [pot] [fillout_factor]
    <!-- http://phoebe-project.org/docs/2.4/tutorials/requiv_crit_contact -->
Surface Gravities
Eccentricity & Volume Conservation
Spots
Eclipse Detection
Passband Luminosity
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