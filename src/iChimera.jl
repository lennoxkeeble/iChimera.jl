module iChimera
include("main.jl")

##### REST OF THE CODE DEFINES AN iChimera TYPE AND THEN RUNS A SHORT INSPIRAL FOR PRECOMPILATION #####
using PrecompileTools: @setup_workload, @compile_workload    # this is a small dependency
using StaticArrays
using ..InclinationMappings
using ..ConstantsOfMotion
using ..Kerr
using ..Inspiral
using ..MultipoleDerivs

const c::Float64 = 2.99792458 * 1e8 # speed of light
const Grav_Newton::Float64 = 6.67430 * 1e-11 # newton's gravitational constant
const Msol::Float64 = (1.988) * 1e30 # one solar mass in kg
const year::Float64 = 365.0 * 24.0 * 60.0 * 60.0 # seconds in a year
const pc::Float64 = 3.08568025e16; # parsec in meters
const kg_to_meters::Float64 = Grav_Newton / c / c;
const obs_distance::Float64 = 1.0;

mutable struct EMRI
    a::Float64 # spin of the (massive) black hole
    p::Float64 # semi-latus rectum
    e::Float64 # eccentricity
    inclination::Float64 # inclination angle in DEGREES
    inclination_type::String # type of inclination. "iota": cos(ι) = Lz / sqrt(Lz^2 + C) (Eq. 25 of arXiv:1109.0572v2) or "theta_inc": θ_inc = π/2 - sign(Lz) * θmin (Eq. 1.2 of arXiv:2401.09577v2)
    θmin::Float64 # minimum polar angle (in radians)
    sign_Lz::Int64 # sign of Lz. 1 for prograde, -1 for retrograde
    mass_ratio::Float64 # mass ratio of the small compact object to the central massive black hole
    lmax_mass_fluxes::Int64 # maximum mass-type multipole moment l mode to include in the flux computation with 2 ≤ lmax ≤ 4
    lmax_current_fluxes::Int64 # maximum current-type multipole moment l mode to include in the flux computation with 1 ≤ lmax ≤ 2 (lmax = 1 excludes any current-type moment and only up to l=2 included at this time)
    lmax_mass_waveform::Int64 # maximum mass-type multipole moment l mode to include in the waveform computation with 2 ≤ lmax ≤ 4
    lmax_current_waveform::Int64 # maximum current-type multipole moment l mode to include in the waveform computation with 1 ≤ lmax ≤ 3 (lmax = 1 excludes any current-type moment and only up to l=3 included at this time)
    coordinates::String # coordinate system for the self-force computation: "harmonic" or "cartesian"
    OnePN::Float64 # 1PN correction to the mass quadrupole moment
    TwoPN::Float64 # 2PN correction to the mass quadrupole moment
    TwoPointFivePN::Float64 # 2.5PN correction to the mass quadrupole moment
    psi0::Float64 # (initial condition) intial radial angle variable
    chi0::Float64 # (initial condition) initial polar angle variable
    phi0::Float64 # (initial condition) initial azimuthal angle variable
    frame::String # waveform frame: "SSB" for solar system barycenter or "Source" for source frame
    ThetaS::Float64 # (waveform — SSB frame) EMRI system polar orientation in solar system barycenter (SSB) frame
    PhiS::Float64 # (waveform — SSB frame) EMRI system azimuthal orientation in SSB frame
    ThetaK::Float64 # (waveform — SSB frame) MBH spin polar orientation in SSB frame
    PhiK::Float64 # (waveform — SSB frame) MBH spin azimuthal orientation in SSB frame
    ThetaObs::Float64 # (waveform — source frame) observer polar orientation
    PhiObs::Float64 # (waveform — source frame) observer azimuthal orientation
    dt_save::Float64 # time interval in seconds between saving data points (e.g., waveform, trajectory, etc.)
    path::String # path for saving output files
    T_secs::Float64 # maximum orbit evolution time in seconds
    M::Float64 # mass of the massive black hole
    reltol::Float64 # relative tolerance for the geodesic ODE solver
    abstol::Float64 # absolute tolerance for the geodesic ODE solver
    compute_SF_frac::Float64 # determines how often the self-force is to be computed as a fraction of the maximum time period: Δt_SF = compute_SF_frac * (2π / min(ω)), where ω are the fundamental frequencies
    save_every::Int64 # save solution to file after every save_every steps
    save_traj::Bool # save BL trajectory
    save_constants::Bool # save constants of motion
    save_fluxes::Bool # save fluxes
    save_gamma::Bool # save gamma factor
    # checks to make sure that the parameters are valid
    EMRI(a, p, e, inclination, inclination_type, θmin, sign_Lz, mass_ratio, lmax_mass_fluxes, lmax_current_fluxes, lmax_mass_waveform, lmax_current_waveform, coordinates, OnePN, TwoPN, TwoPointFivePN, psi0, chi0, phi0, frame, ThetaS,
        PhiS, ThetaK, PhiK, ThetaObs, PhiObs, dt_save, path, T_secs, M, reltol, abstol, compute_SF_frac, save_every, save_traj, save_constants, save_fluxes, save_gamma) = begin
    if a < 0.0 || a > 1.0
        error("Spin parameter 'a' must be between 0 and 1.")
    elseif p <= 0.0
        error("Semi-latus rectum 'p' must be positive.")
    elseif e < 0.0 || e >= 1.0
        error("Eccentricity 'e' must be between 0 and 1.")
    elseif inclination < 0.0 || inclination > 180.0
        error("Inclination 'inclination' must be between 0 and 180 degrees.")
    elseif inclination_type != "iota" && inclination_type != "theta_inc"
        error("Inclination type 'inclination_type' must be either 'iota' or 'theta_inc'.")
    elseif θmin < 0.0 || θmin > π
        error("Minimum polar angle 'θmin' must be between 0 and π radians.")
    elseif sign_Lz != 1 && sign_Lz != -1
        error("Sign of Lz 'sign_Lz' must be either 1 (prograde) or -1 (retrograde).")
    elseif mass_ratio <= 0.0
        error("Mass ratio 'mass_ratio' must be positive.")
    elseif lmax_mass_fluxes < 2 || lmax_mass_fluxes > 3
        error("Maximum mode 'lmax_mass_fluxes' must be between 2 and 3.")
    elseif lmax_current_fluxes < 1 || lmax_current_fluxes > 2
        error("Maximum mode 'lmax_current_fluxes' must be between 1 and 2.")
    elseif lmax_mass_waveform < 2 || lmax_mass_waveform > 4
        error("Maximum mode 'lmax_mass_waveform' must be between 2 and 4.")
    elseif lmax_current_waveform < 1 || lmax_current_waveform > 3
        error("Maximum mode 'lmax_current_waveform' must be between 1 and 3.")
    elseif coordinates != "harmonic" && coordinates != "cartesian"
        error("Coordinate system must be either 'harmonic' or 'cartesian'.")
    elseif OnePN != 0.0 && OnePN != 1.0
        error("1PN correction 'OnePN' must be either 0.0 (off) or 1.0 (on).")
    elseif TwoPN != 0.0 && TwoPN != 1.0
        error("2PN correction 'TwoPN' must be either 0.0 (off) or 1.0 (on).")
    elseif TwoPointFivePN != 0.0 && TwoPointFivePN != 1.0
        error("2.5PN correction 'TwoPointFivePN' must be either 0.0 (off) or 1.0 (on).")
    elseif frame != "SSB" && frame != "Source"
        error("Waveform frame 'frame' must be either 'SSB' or 'Source'.")
    elseif ThetaS < 0.0 || ThetaS > 180.0
        error("ThetaS 'ThetaS' must be between 0 and 180 degrees.")
    elseif PhiS < 0.0 || PhiS > 360.0
        error("PhiS 'PhiS' must be between 0 and 360 degrees.")
    elseif ThetaK < 0.0 || ThetaK > 180.0
        error("ThetaK 'ThetaK' must be between 0 and 180 degrees.")
    elseif PhiK < 0.0 || PhiK > 360.0
        error("PhiK 'PhiK' must be between 0 and 360 degrees.")
    elseif ThetaObs < 0.0 || ThetaObs > 180.0
        error("ThetaObs 'ThetaObs' must be between 0 and 180 degrees.")
    elseif PhiObs < 0.0 || PhiObs > 360.0
        error("PhiObs 'PhiObs' must be between 0 and 360 degrees.")
    elseif dt_save <= 0.0
        error("Time interval between saving data points 'dt_save' must be positive.")
    elseif T_secs <= 0.0
        error("Maximum orbit evolution time 'T_secs' must be positive.")
    elseif M <= 0.0
        error("Mass of the massive black hole 'M' must be positive.")
    elseif reltol <= 0.0
        error("Relative tolerance 'reltol' must be positive.")
    elseif abstol <= 0.0
        error("Absolute tolerance 'abstol' must be positive.")
    elseif compute_SF_frac <= 0.0
        error("Compute SF fraction 'compute_SF_frac' must be positive.")
    elseif save_every <= 0
        error("Save every 'save_every' must be positive.")
    else
        new(a, p, e, inclination, inclination_type, θmin, 
        sign_Lz, mass_ratio, lmax_mass_fluxes, lmax_current_fluxes, lmax_mass_waveform, lmax_current_waveform, coordinates, OnePN, TwoPN, TwoPointFivePN, psi0, chi0, phi0, frame, ThetaS, PhiS, ThetaK, PhiK, ThetaObs, PhiObs, dt_save, path, T_secs, M, reltol, abstol, compute_SF_frac, save_every, save_traj, save_constants, save_fluxes, save_gamma)
    end
    end
end

function compute_theta_min(a::Float64, p::Float64, e::Float64, inclination::Float64, inclination_type::String, sign_Lz::Int64)
    if inclination == 0.0
        return π/2
    elseif inclination_type == "iota"
        return InclinationMappings.iota_to_theta_min(a, p, e, inclination)
    elseif inclination_type == "theta_inc"
        return InclinationMappings.theta_inc_to_theta_min(inclination, sign_Lz)
    else
        error("Invalid inclination type. Use 'iota' or 'theta_inc'.")
    end
end

function EMRI(
    a::Float64,
    p::Float64,
    e::Float64,
    inclination::Float64,
    inclination_type::String,
    sign_Lz::Int64,
    mass_ratio::Float64,
    lmax_mass_fluxes::Int64,
    lmax_current_fluxes::Int64,
    lmax_mass_waveform::Int64,
    lmax_current_waveform::Int64,
    coordinates::String,
    OnePN::Float64,
    TwoPN::Float64,
    TwoPointFivePN::Float64,
    psi0::Float64,
    chi0::Float64,
    phi0::Float64,
    frame::String,
    ThetaS::Float64,
    PhiS::Float64,
    ThetaK::Float64,
    PhiK::Float64,
    ThetaObs::Float64,
    PhiObs::Float64,
    dt_save::Float64,
    path::String,
    T_secs::Float64,
    M::Float64,
    reltol::Float64,
    abstol::Float64,
    compute_SF_frac::Float64,
    save_every::Int64,
    save_traj::Bool,
    save_constants::Bool,
    save_fluxes::Bool,
    save_gamma::Bool)
    # checks to make sure that the parameters are valid
    if a < 0.0 || a > 1.0
        error("Spin parameter 'a' must be between 0 and 1.")
    elseif p <= 0.0
        error("Semi-latus rectum 'p' must be positive.")
    elseif e < 0.0 || e >= 1.0
        error("Eccentricity 'e' must be between 0 and 1.")
    elseif inclination < 0.0 || inclination > 180.0
        error("Inclination 'inclination' must be between 0 and 180 degrees.")
    elseif inclination_type != "iota" && inclination_type != "theta_inc"
        error("Inclination type 'inclination_type' must be either 'iota' or 'theta_inc'.")
    elseif sign_Lz != 1 && sign_Lz != -1
        error("Sign of Lz 'sign_Lz' must be either 1 (prograde) or -1 (retrograde).")
    elseif mass_ratio <= 0.0
        error("Mass ratio 'mass_ratio' must be positive.")
    elseif lmax_mass_fluxes < 2 || lmax_mass_fluxes > 3
        error("Maximum mode 'lmax_mass_fluxes' must be between 2 and 3.")
    elseif lmax_current_fluxes < 1 || lmax_current_fluxes > 2
        error("Maximum mode 'lmax_current_fluxes' must be between 1 and 2.")
    elseif lmax_mass_waveform < 2 || lmax_mass_waveform > 4
        error("Maximum mode 'lmax_mass_waveform' must be between 2 and 4.")
    elseif lmax_current_waveform < 1 || lmax_current_waveform > 3
        error("Maximum mode 'lmax_current_waveform' must be between 1 and 3.")
    elseif coordinates != "harmonic" && coordinates != "cartesian"
        error("Coordinate system must be either 'harmonic' or 'cartesian'.")
    elseif OnePN != 0.0 && OnePN != 1.0
        error("1PN correction 'OnePN' must be either 0.0 (off) or 1.0 (on).")
    elseif TwoPN != 0.0 && TwoPN != 1.0
        error("2PN correction 'TwoPN' must be either 0.0 (off) or 1.0 (on).")
    elseif TwoPointFivePN != 0.0 && TwoPointFivePN != 1.0
        error("2.5PN correction 'TwoPointFivePN' must be either 0.0 (off) or 1.0 (on).")
    elseif frame != "SSB" && frame != "Source"
        error("Waveform frame 'frame' must be either 'SSB' or 'Source'.")
    elseif ThetaS < 0.0 || ThetaS > 180.0
        error("ThetaS 'ThetaS' must be between 0 and 180 degrees.")
    elseif PhiS < 0.0 || PhiS > 360.0
        error("PhiS 'PhiS' must be between 0 and 360 degrees.")
    elseif ThetaK < 0.0 || ThetaK > 180.0
        error("ThetaK 'ThetaK' must be between 0 and 180 degrees.")
    elseif PhiK < 0.0 || PhiK > 360.0
        error("PhiK 'PhiK' must be between 0 and 360 degrees.")
    elseif ThetaObs < 0.0 || ThetaObs > 180.0
        error("ThetaObs 'ThetaObs' must be between 0 and 180 degrees.")
    elseif PhiObs < 0.0 || PhiObs > 360.0
        error("PhiObs 'PhiObs' must be between 0 and 360 degrees.")
    elseif dt_save <= 0.0
        error("Time interval between saving data points 'dt_save' must be positive.")
    elseif T_secs <= 0.0
        error("Maximum orbit evolution time 'T_secs' must be positive.")
    elseif M <= 0.0
        error("Mass of the massive black hole 'M' must be positive.")
    elseif reltol <= 0.0
        error("Relative tolerance 'reltol' must be positive.")
    elseif abstol <= 0.0
        error("Absolute tolerance 'abstol' must be positive.")
    elseif compute_SF_frac <= 0.0
        error("Compute SF fraction 'compute_SF_frac' must be positive.")
    elseif save_every <= 0
        error("Save every 'save_every' must be positive.")
    end
    theta_min = compute_theta_min(a, p, e, inclination, inclination_type, sign_Lz)
    return EMRI(a, p, e, inclination, inclination_type, theta_min, sign_Lz, mass_ratio, lmax_mass_fluxes, lmax_current_fluxes, lmax_mass_waveform, lmax_current_waveform, coordinates, OnePN, TwoPN, TwoPointFivePN, psi0, chi0, phi0, frame, ThetaS,
        PhiS, ThetaK, PhiK, ThetaObs, PhiObs, dt_save, path, T_secs, M, reltol, abstol, compute_SF_frac, save_every, save_traj, save_constants, save_fluxes, save_gamma)
end

function compute_inspiral(emri; JIT::Bool = false, rr_model::Union{Symbol, AbstractString}=:chimera)
    a = emri.a
    p = emri.p
    e = emri.e
    θmin = emri.θmin
    sign_Lz = emri.sign_Lz
    inclination = emri.inclination
    Mass_MBH = emri.M
    t_max_secs = emri.T_secs
    compute_SF_frac = emri.compute_SF_frac


    EE, LL, QQ, CC = ConstantsOfMotion.compute_ELC(a, p, e, θmin, sign_Lz);
    rplus = Kerr.KerrMetric.rplus(a); rminus = Kerr.KerrMetric.rminus(a);

    # Mino time frequencies
    ω = ConstantsOfMotion.KerrFreqs(a, p, e, θmin, EE, LL, QQ, CC, rplus, rminus);

    # BL time frequencies
    Ω = ω[1:3]/ω[4]; Ωr, Ωθ, Ωϕ = Ω;

    ### evolution time ###
    MtoSecs = Mass_MBH * Grav_Newton / c^3; # conversion from t(M) -> t(s)
    t_max_M = t_max_secs / MtoSecs; # units of M
    dt_save_M = emri.dt_save / MtoSecs; # units of M

    
    if e != 0.0 && inclination != 0.0
        compute_fluxes = compute_SF_frac * minimum(@. 2π /Ω[1:3])
    # eccentric equatorial
    elseif e != 0.0
        compute_fluxes = compute_SF_frac * minimum(@. 2π /[Ω[1], Ω[3]])
    # circular inclined
    elseif inclination != 0.0
        compute_fluxes = compute_SF_frac * minimum(@. 2π /Ω[2:3])
    end

    Inspiral.compute_inspiral(emri.a, emri.p, emri.e, emri.θmin, emri.sign_Lz, emri.mass_ratio, emri.psi0, emri.chi0, emri.phi0, compute_fluxes, t_max_M, dt_save_M, emri.save_every, emri.reltol, emri.abstol, emri.OnePN, emri.TwoPN, emri.TwoPointFivePN, emri.coordinates; data_path=emri.path, JIT=JIT, lmax_mass_fluxes=emri.lmax_mass_fluxes, lmax_current_fluxes=emri.lmax_current_fluxes, save_traj=emri.save_traj, save_constants=emri.save_constants, save_fluxes=emri.save_fluxes, save_gamma=emri.save_gamma, rr_model=rr_model)
end

function compute_waveform(emri; rr_model::Union{Symbol, AbstractString}=:chimera)
    if emri.frame == "SSB"
        Inspiral.compute_waveform(obs_distance, emri.ThetaS, emri.PhiS, emri.ThetaK, emri.PhiK, emri.a, emri.p, emri.e, emri.θmin, emri.mass_ratio, emri.psi0, emri.chi0, emri.phi0, emri.lmax_mass_fluxes, emri.lmax_current_fluxes, emri.lmax_mass_waveform, emri.lmax_current_waveform, emri.coordinates, emri.path; rr_model=rr_model);
    elseif emri.frame == "Source"
        Inspiral.compute_waveform(obs_distance, emri.ThetaObs, emri.PhiObs, emri.a, emri.p, emri.e, emri.θmin, emri.mass_ratio, emri.psi0, emri.chi0, emri.phi0, emri.lmax_mass_fluxes, emri.lmax_current_fluxes, emri.lmax_mass_waveform, emri.lmax_current_waveform, emri.coordinates, emri.path; rr_model=rr_model);
    end
end

function load_trajectory(emri; rr_model::Union{Symbol, AbstractString}=:chimera)
    fname = Inspiral.solution_fname(emri.a, emri.p, emri.e, emri.θmin, emri.mass_ratio, emri.psi0, emri.chi0, emri.phi0, emri.lmax_mass_fluxes, emri.lmax_current_fluxes, emri.coordinates, emri.path; rr_model=rr_model)
    return Inspiral.load_trajectory(fname)
end

function load_constants_of_motion(emri; rr_model::Union{Symbol, AbstractString}=:chimera)
    fname = Inspiral.solution_fname(emri.a, emri.p, emri.e, emri.θmin, emri.mass_ratio, emri.psi0, emri.chi0, emri.phi0, emri.lmax_mass_fluxes, emri.lmax_current_fluxes, emri.coordinates, emri.path; rr_model=rr_model)
    return Inspiral.load_constants_of_motion(fname)
end

function load_fluxes(emri; rr_model::Union{Symbol, AbstractString}=:chimera)
    fname = Inspiral.solution_fname(emri.a, emri.p, emri.e, emri.θmin, emri.mass_ratio, emri.psi0, emri.chi0, emri.phi0, emri.lmax_mass_fluxes, emri.lmax_current_fluxes, emri.coordinates, emri.path; rr_model=rr_model)
    return Inspiral.load_fluxes(fname)
end

function load_waveform(emri; rr_model::Union{Symbol, AbstractString}=:chimera)
    if emri.frame == "SSB"
        return Inspiral.load_waveform(obs_distance, emri.ThetaS, emri.PhiS, emri.ThetaK, emri.PhiK, emri.a, emri.p, emri.e, emri.θmin, emri.mass_ratio, emri.psi0, emri.chi0, emri.phi0, emri.lmax_mass_fluxes, emri.lmax_current_fluxes, emri.lmax_mass_waveform, emri.lmax_current_waveform, emri.coordinates, emri.path; rr_model=rr_model);
    elseif emri.frame == "Source"
        return Inspiral.load_waveform(obs_distance, emri.ThetaObs, emri.PhiObs, emri.a, emri.p, emri.e, emri.θmin, emri.mass_ratio, emri.psi0, emri.chi0, emri.phi0, emri.lmax_mass_fluxes, emri.lmax_current_fluxes, emri.lmax_mass_waveform, emri.lmax_current_waveform, emri.coordinates, emri.path; rr_model=rr_model);
    end
end

# precompilation
@setup_workload begin
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

    OnePN = 1.0
    TwoPN = 1.0
    TwoPointFivePN = 1.0

    coordinates = "cartesian";

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
    reltol =  1e-14; # relative tolerance for the geodesic solver
    abstol =  1e-14; # absolute tolerance for the geodesic solver

    # precompile mass quadrupole calls

    @compile_workload begin
        emri = iChimera.EMRI(a, p, e, inclination, inclination_type, sign_Lz, mass_ratio, lmax_mass_fluxes, lmax_current_fluxes, lmax_mass_waveform, lmax_current_waveform, coordinates, OnePN, TwoPN, TwoPointFivePN, psi0, chi0, phi0, frame, ThetaS, PhiS, ThetaK, PhiK, ThetaObs, PhiObs, dt_save, data_path, t_max_secs, Mass_MBH, reltol, abstol, compute_SF_frac, save_every, save_traj, save_constants, save_fluxes, save_gamma);
        @time iChimera.compute_inspiral(emri; JIT = true);
    end
end

end
