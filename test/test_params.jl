using iChimera

# (initial) orbital parameters
a = 0.98; # spin of the (massive) black hole
e = 0.6; # eccentricity
mass_ratio = 1e-5; # mass ratio of the small compact object to the central massive black hole
p = 7.0; # semi-latus rectum

inclination = 57.39; # inclination angle (in degrees)
sign_Lz = inclination < 90.0 ? 1 : -1; # sign of z-component of angular momentum: +1 for prograde, -1 for retrograde
inclination_type = "iota"; # type of inclination. "iota": cos(ι) = Lz / sqrt(Lz^2 + C) (Eq. 25 of arXiv:1109.0572v2) or "theta_inc": θ_inc = π/2 - sign(Lz) * θmin (Eq. 1.2 of arXiv:2401.09577v2)

t_max_secs = (10^-3) * iChimera.year / 3.; # maximum orbit evolution time in seconds
Mass_MBH = 1e6 * iChimera.Msol; # mass of the massive black hole

dt_save = 5.0; # time interval in seconds between saving data points (e.g., waveform, trajectory, etc.)
save_every = 1000; # save solution to file after every save_every steps

# initial angle variables
psi0 = 0.1; # (initial condition) intial radial angle variable
chi0 = 0.2; # (initial condition) initial polar angle variable
phi0 = 0.3; # (initial condition) initial azimuthal angle variable

# waveform parameters
obs_distance_Gpc = 0.5;  # radial distance to observer in Gpc
frame = "SSB"; # frame in which the waveform is computed: "SSB" for solar system barycenter or "Source" for source frame
ThetaS = 10.0; # (waveform — SSB) EMRI system polar orientation in solar system barycenter (SSB) frame (degrees)
PhiS = 5.0; # (waveform — SSB) EMRI system azimuthal orientation in SSB frame (degrees)
ThetaK = 6.0; # (waveform — SSB) MBH spin polar orientation in SSB frame (degrees)
PhiK = 8.0; # (waveform — SSB) MBH spin azimuthal orientation in SSB frame (degrees)
ThetaObs = 50.0; # (waveform — source) observer polar orientation in source frame (degrees)
PhiObs = 20.0; # (waveform — source) observer azimuthal orientation in source frame (degrees)

# chose multipole moments to include in flux and waveform
lmax_mass_fluxes = 3 # maximum mass-type multipole moment l mode to include in the flux computation with 2 ≤ lmax ≤ 3 (mass quadrupole and octupole)
lmax_current_fluxes = 2 # maximum current-type multipole moment l mode to include in the flux computation with 1 ≤ lmax ≤ 2 (lmax = 1 excludes any current-type moment and only up to current quadrupole included at this time)

lmax_mass_waveform = 4 # maximum mass-type multipole moment l mode to include in the waveform computation with 2 ≤ lmax ≤ 4 (mass quadrupole, octupole, and hexadecapole)
lmax_current_waveform = 3 # maximum current-type multipole moment l mode to include in the waveform computation with 1 ≤ lmax ≤ 3 (lmax = 1 excludes any current-type moment and only current quadrupole and octupole included at this time)

# data saving options
save_traj = true; # save physical BL trajectory
save_SF = false; # save BL self-force
save_constants = true; # save constants of motion
save_fluxes = true; # save fluxes
save_gamma = true; # save gamma factor

# file paths
results_path = "./Results";
data_path=results_path * "/Data/"; # path for saving output files
mkpath(data_path);

# compute multiplicative factor to convert GW strain into SI units
MtoSecs = Mass_MBH * iChimera.Grav_Newton / iChimera.c^3
M_to_meters = iChimera.kg_to_meters * Mass_MBH;
meters_to_M = 1.0 / M_to_meters;
Gpc_to_M = 1.0e9 * iChimera.pc * meters_to_M; # 1 Gpc in units of M
observer_distance_M = obs_distance_Gpc * Gpc_to_M;
strain_to_SI = mass_ratio / observer_distance_M;

compute_SF_frac = 0.01 # determines how often the self-force is to be computed as a fraction of the maximum time period: Δt_SF = compute_SF_frac * (2π / min(ω)), where ω are the fundamental frequencies

# geodesic ODE solver options
reltol =  1e-16; # relative tolerance for the geodesic solver
abstol =  1e-16; # absolute tolerance for the geodesic solver

OnePN = 1.0
TwoPN = 1.0
TwoPointFivePN = 1.0

coordinates = "cartesian";

# intialize EMRI type
emri = iChimera.EMRI(a, p, e, inclination, inclination_type, sign_Lz, mass_ratio, lmax_mass_fluxes, lmax_current_fluxes, lmax_mass_waveform, lmax_current_waveform, coordinates, OnePN, TwoPN, TwoPointFivePN, psi0, chi0, phi0, frame, ThetaS, PhiS, ThetaK, PhiK, ThetaObs, PhiObs, dt_save, data_path, t_max_secs, Mass_MBH, reltol, abstol, compute_SF_frac, save_every, save_traj, save_constants, save_fluxes, save_gamma);