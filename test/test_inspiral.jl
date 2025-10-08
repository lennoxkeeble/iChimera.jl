include("../src/main.jl")
using iChimera
using HDF5
using LinearAlgebra
using Combinatorics
using StaticArrays
using Printf
using DifferentialEquations
import ..HDF5Helper: create_file_group!, create_dataset!, append_data!
using ..SymmetricTensors
using ..Kerr
using ..ConstantsOfMotion
using ..BLTimeGeodesics
using ..SelfAccelerationHarmonic
using ..EvolveConstants
using ..Waveform
using ..HarmonicCoordDerivs
using ..CoordinateDerivs
using ..MultipoleDerivs

Z_1(a) = 1 + (1 - a^2 / 1.0)^(1/3) * ((1 + a)^(1/3) + (1 - a)^(1/3))
Z_2(a) = sqrt(3 * a^2 + Z_1(a)^2)
LSO_r(a) = (3 + Z_2(a) - sqrt((3 - Z_1(a)) * (3 + Z_1(a) * 2 * Z_2(a))))   # retrograde LSO
LSO_p(a) = (3 + Z_2(a) + sqrt((3 - Z_1(a)) * (3 + Z_1(a) * 2 * Z_2(a))))   # prograde LSO


a = 0.8;
p = 10.0;
e = 0.5;
θmin = π/3;
sign_Lz = 1;
q = 1e-5;
psi_0 = 0.0;
chi_0 = 0.0;
phi_0 = 0.0;
compute_SF = 0.1;
tInspiral = 10.0;
dt_save = 0.5;
save_every = 100;
reltol = 1e-15;
abstol = 1e-15;
data_path = "Data/"
mkpath(data_path)
JIT = false
lmax_mass_fluxes = 3
lmax_current_fluxes = 2
save_traj = false
save_constants = false
save_fluxes = false
save_gamma = false
maxiters = Int(1e8)


# create solution file
sol_filename=iChimera.Inspiral.solution_fname(a, p, e, θmin, q, psi_0, chi_0, phi_0, lmax_mass_fluxes, lmax_current_fluxes, data_path)

if isfile(sol_filename)
    rm(sol_filename)
end

if JIT
    tInspiral = 20.0 # dummy run for Δt = 20M
end

file = h5open(sol_filename, "w")

# second argument is chunk_size. Since each successive geodesic piece overlap at the end of the first and bgeinning of the second, we must manually save this point only once to avoid repeats in the data 
Inspiral.initialize_solution_file!(file, save_every, save_traj, save_constants, save_fluxes, save_gamma)

# initialize data arrays for trajectory and multipole moments which will be used for post-processing waveform computation
idx_save_1 = 1;
time = zeros(save_every)
r = zeros(save_every)
theta = zeros(save_every)
phi = zeros(save_every)
gamma = zeros(save_every)
Mij2_data = [zeros(save_every) for i=1:3, j=1:3]
Sij2_data = [zeros(save_every) for i=1:3, j=1:3];
Mijk3_data = [zeros(save_every) for i=1:3, j=1:3, k=1:3];
Sijk3_data = [zeros(save_every) for i=1:3, j=1:3, k=1:3];
Mijkl4_data = [zeros(save_every) for i=1:3, j=1:3, k=1:3, l=1:3];

# create arrays to store multipole moments necessary for waveform computation computed in a given time step (the above are used to store over multiple time steps)
Mij2 = @MArray zeros(3, 3);
Mijk3 = @MArray zeros(3, 3, 3);
Mijkl4 = @MArray zeros(3, 3, 3, 3);
Sij2 = @MArray zeros(3, 3);
Sijk3 = @MArray zeros(3, 3, 3);

# initialize derivative arrays
xBL = @MArray zeros(3); vBL = @MArray zeros(3); aBL = @MArray zeros(3);

dxBL_dt = @MArray zeros(3); d2xBL_dt = @MArray zeros(3); d3xBL_dt = @MArray zeros(3); d4xBL_dt = @MArray zeros(3);
d5xBL_dt = @MArray zeros(3); d6xBL_dt = @MArray zeros(3); d7xBL_dt = @MArray zeros(3); d8xBL_dt = @MArray zeros(3);

dx_dλ = @MArray zeros(3); d2x_dλ = @MArray zeros(3); d3x_dλ = @MArray zeros(3); d4x_dλ = @MArray zeros(3);
d5x_dλ = @MArray zeros(3); d6x_dλ = @MArray zeros(3); d7x_dλ = @MArray zeros(3); d8x_dλ = @MArray zeros(3);

xH = @MArray zeros(3); dxH_dt = @MArray zeros(3); d2xH_dt = @MArray zeros(3); d3xH_dt = @MArray zeros(3); d4xH_dt = @MArray zeros(3);
d5xH_dt = @MArray zeros(3); d6xH_dt = @MArray zeros(3); d7xH_dt = @MArray zeros(3); d8xH_dt = @MArray zeros(3);

vH = @MArray zeros(3);
aH = @MArray zeros(3);

# arrays for self-force computation
aSF_BL = @MArray zeros(4)
aSF_H = @MArray zeros(4)
Virr = @MArray zeros(3)
∂Vrr_∂a = @MArray zeros(3)
∂Virr_∂t = @MArray zeros(3)
∂Virr_∂a = @MArray zeros(3, 3)

# compute apastron
ra = p / (1 - e);

# calculate integrals of motion from orbital parameters
EEi, LLi, QQi, CCi = ConstantsOfMotion.compute_ELC(a, p, e, θmin, sign_Lz)   

# store orbital params in arrays
idx_save_2 = 1
t_Fluxes = zeros(save_every);
E_arr = zeros(save_every);
E_dot_arr = zeros(save_every);
L_arr = zeros(save_every); 
L_dot_arr = zeros(save_every);
C_arr = zeros(save_every);
C_dot_arr = zeros(save_every);
Q_arr = zeros(save_every);
Q_dot_arr = zeros(save_every);
p_arr = zeros(save_every);
e_arr = zeros(save_every);
θmin_arr = zeros(save_every);

E_t = EEi; 
dE_dt = 0.;
L_t = LLi; 
dL_dt = 0.;
C_t = CCi;
dC_dt = 0.;
Q_t = QQi
dQ_dt = 0.;
p_t = p;
e_t = e;
θmin_t = θmin;

rplus = Kerr.KerrMetric.rplus(a); rminus = Kerr.KerrMetric.rminus(a);

# initial condition for Kerr geodesic trajectory
t0 = 0.0
rLSO = Inspiral.LSO_p(a)

# initialize ODE problem
E, L, Q, C, ra, p3, p4, zp, zm = BLTimeGeodesics.compute_ODE_params(a, p, e, θmin, sign_Lz);

params = @SArray [a, E, L, p, e, θmin, p3, p4, zp, zm];
ics = @SArray[psi_0, chi_0, phi_0];

# initial conditions for Kerr geodesic trajectory
tspan = (0.0, tInspiral);

prob = e == 0.0 ? ODEProblem(BLTimeGeodesics.HJ_Eqns_circular, ics, tspan, params) : ODEProblem(BLTimeGeodesics.HJ_Eqns, ics, tspan, params);

# times at which the integrator should be stopped (consists of times at which the self-force must be computed and times at which the waveform must be computed)
times_SF = range(start = compute_SF, stop = tInspiral, step = compute_SF) |> collect;
times_WF = range(start = dt_save, stop = tInspiral, step = dt_save) |> collect;

# initialize integrator
Δti = compute_SF / 1000;
integrator = init(prob, AutoTsit5(RK4()), adaptive=true, dt = Δti, reltol = reltol, abstol = abstol, maxiters = maxiters)

#### SAVE SOLUTION AT TIME t = 0 ####;
## COMPUTE HIGH-ORDER TIME DERIVATIVES OF THE TRAJECTORY ##
# EXTRACT GEODESIC TRAJECTORY FROM INTEGRATOR
tt, rr, θθ, ϕϕ, r_dot, θ_dot, ϕ_dot, r_ddot, θ_ddot, ϕ_ddot, dt_dτ, psi, chi = Inspiral.compute_geodesic_arrays(integrator, a, E_t, L_t, p_t, e_t, θmin_t, p3, p4, zp, zm)

# COMPUTE BL COORDINATE DERIVATIVES
xBL[1] = rr; xBL[2] = θθ; xBL[3] = ϕϕ;
vBL[1] = r_dot; vBL[2] = θ_dot; vBL[3] = ϕ_dot;
aBL[1] = r_ddot; aBL[2] = θ_ddot; aBL[3] = ϕ_ddot;
CoordinateDerivs.ComputeDerivs!(xBL, sign(vBL[1]), sign(vBL[2]), dxBL_dt, d2xBL_dt, d3xBL_dt, d4xBL_dt, d5xBL_dt, d6xBL_dt, d7xBL_dt, d8xBL_dt, dx_dλ, d2x_dλ, d3x_dλ, d4x_dλ, d5x_dλ, d6x_dλ, d7x_dλ, d8x_dλ, a, E_t, L_t, C_t);

# COMPUTE HARMONIC COORDINATE DERIVATIVES
HarmonicCoordDerivs.compute_harmonic_derivs!(xBL, dxBL_dt, d2xBL_dt, d3xBL_dt, d4xBL_dt, d5xBL_dt, d6xBL_dt, d7xBL_dt, d8xBL_dt, xH, dxH_dt, d2xH_dt, d3xH_dt, d4xH_dt, d5xH_dt, d6xH_dt, d7xH_dt, d8xH_dt, a);

# UPDATE MULTIPOLE MOMENTS AND TRAJECTORY ARRAYS
MultipoleDerivs.compute_WF_moments!(q, Mij2, Mijk3, Mijkl4, Sij2, Sijk3, xH, dxH_dt, d2xH_dt, d3xH_dt, d4xH_dt); # note that all implemented moments are calculated for the waveform. In the post-processing step, one can choose which moments to include
Inspiral.update_waveform_arrays!(idx_save_1, Mij2_data, Sij2_data, Mijk3_data, Sijk3_data, Mijkl4_data, Mij2, Sij2, Mijk3, Sijk3, Mijkl4)
Inspiral.update_trajectory_arrays!(integrator, idx_save_1, time, r, theta, phi, gamma, a, E_t, L_t, p_t, e_t, θmin_t, p3, p4, zp, zm, save_traj, save_gamma)
idx_save_1 += 1

# also save constants of motion
t_Fluxes[idx_save_2] = t0;
E_arr[idx_save_2] = E_t;
E_dot_arr[idx_save_2] = dE_dt;
L_arr[idx_save_2] = L_t; 
L_dot_arr[idx_save_2] = dL_dt;
C_arr[idx_save_2] = C_t;
C_dot_arr[idx_save_2] = dC_dt;
Q_arr[idx_save_2] = Q_t;
Q_dot_arr[idx_save_2] = dQ_dt;
p_arr[idx_save_2] = p_t;
e_arr[idx_save_2] = e_t;
θmin_arr[idx_save_2] = θmin_t;
idx_save_2 += 1

WF_step = false
SF_step = true

t0 = 0.0
# while integrator.t < tInspiral
    # print("Completion: $(round(100 * t0/tInspiral; digits=5))%   \r")
    # flush(stdout)
    
    # # successively step the integrator to the next time in either times_SF or times_WF, whichever is smaller
    # if length(times_SF) == 0 && length(times_WF) == 0
    #     break
    # elseif length(times_WF) == 0
    #     tF = times_SF[1]
    #     popfirst!(times_SF)
    #     SF_step = true
    #     WF_step = false
    # elseif length(times_SF) == 0
    #     tF = times_WF[1]
    #     popfirst!(times_WF)
    #     WF_step = true
    #     SF_step = false
    # elseif times_SF[1] < times_WF[1]
    #     tF = times_SF[1]
    #     popfirst!(times_SF)
    #     SF_step = true
    #     WF_step = false
    # else
    #     tF = times_WF[1]
    #     popfirst!(times_WF)
    #     WF_step = true
    #     SF_step = false
    # end

    # time_step = tF - integrator.t
    # step!(integrator, time_step, true)


## COMPUTE HIGH-ORDER TIME DERIVATIVES OF THE TRAJECTORY ##
# EXTRACT GEODESIC TRAJECTORY
tt, rr, θθ, ϕϕ, r_dot, θ_dot, ϕ_dot, r_ddot, θ_ddot, ϕ_ddot, dt_dτ, psi, chi = Inspiral.compute_geodesic_arrays(integrator, a, E_t, L_t, p_t, e_t, θmin_t, p3, p4, zp, zm)

# COMPUTE BL COORDINATE DERIVATIVES
xBL[1] = rr; xBL[2] = θθ; xBL[3] = ϕϕ;
vBL[1] = r_dot; vBL[2] = θ_dot; vBL[3] = ϕ_dot;
aBL[1] = r_ddot; aBL[2] = θ_ddot; aBL[3] = ϕ_ddot;
CoordinateDerivs.ComputeDerivs!(xBL, sign(vBL[1]), sign(vBL[2]), dxBL_dt, d2xBL_dt, d3xBL_dt, d4xBL_dt, d5xBL_dt, d6xBL_dt, d7xBL_dt, d8xBL_dt, dx_dλ, d2x_dλ, d3x_dλ, d4x_dλ, d5x_dλ, d6x_dλ, d7x_dλ, d8x_dλ, a, E_t, L_t, C_t);

# COMPUTE HARMONIC COORDINATE DERIVATIVES
HarmonicCoordDerivs.compute_harmonic_derivs!(xBL, dxBL_dt, d2xBL_dt, d3xBL_dt, d4xBL_dt, d5xBL_dt, d6xBL_dt, d7xBL_dt, d8xBL_dt, xH, dxH_dt, d2xH_dt, d3xH_dt, d4xH_dt, d5xH_dt, d6xH_dt, d7xH_dt, d8xH_dt, a);


## COMPUTE THE FLUXES AND UPDATE ORBITAL CONSTANTS OR COMPUTE WAVEFORM MOMENTS, DEPENDING ON WHETHER IT IS A SELF-FORCE OR WAVEFORM STEP ##
# if WF_step
#     if idx_save_1 == save_every + 1
#         Inspiral.save_traj!(file, save_every, time, r, theta, phi, gamma, save_traj, save_gamma)
#         Inspiral.save_moments!(file, save_every, Mij2_data, Sij2_data, Mijk3_data, Sijk3_data, Mijkl4_data)
#         idx_save_1 = 1
#         flush(file)
#     end

#     MultipoleDerivs.compute_WF_moments!(q, Mij2, Mijk3, Mijkl4, Sij2, Sijk3, xH, dxH_dt, d2xH_dt, d3xH_dt, d4xH_dt); # note that all implemented moments are calculated for the waveform. In the post-processing step, one can choose which moments to include
#     update_waveform_arrays!(idx_save_1, Mij2_data, Sij2_data, Mijk3_data, Sijk3_data, Mijkl4_data, Mij2, Sij2, Mijk3, Sijk3, Mijkl4)
#     update_trajectory_arrays!(integrator, idx_save_1, time, r, theta, phi, gamma, a, E_t, L_t, p_t, e_t, θmin_t, p3, p4, zp, zm, save_traj, save_gamma)
#     idx_save_1 += 1
# else            
##### COMPUTE SELF-FORCE AND SELECT OUT DESIRED MULTIPOLE MOMENTS TO CONTRIBUTE TO SELF-FORCE COMPUTATION (MASS QUADRUPOLE ALWAYS INCLUDED. lmax_mass_fluxes=3 INCLUDES MASS OCTUPOLE, lmax_current_fluxes=2 INCLUDES CURRENT QUADRUPOLE) #####
vH .= dxH_dt
aH .= d2xH_dt
rH = SelfAccelerationHarmonic.norm_3d(xH);
v = SelfAccelerationHarmonic.norm_3d(vH);
SelfAccelerationHarmonic.aRRα(aSF_H, aSF_BL, xH, v, vH, xBL, rH, a, Virr, ∂Vrr_∂a, ∂Virr_∂t, ∂Virr_∂a, xH, dxH_dt, d2xH_dt, d3xH_dt, d4xH_dt, d5xH_dt, d6xH_dt, d7xH_dt, d8xH_dt, q, lmax_mass_fluxes, lmax_current_fluxes)

# test two different flux computations
# dt_flux = tF - t0
dt_flux = 0.01
__, dE_dt, __, dL_dt, __, dQ_dt, __, dC_dt, __, __, __ = EvolveConstants.Evolve_BL(dt_flux, a, rr, θθ, ϕϕ, dt_dτ, r_dot, θ_dot, ϕ_dot, aSF_BL, E_t, dE_dt, L_t, dL_dt, Q_t, dQ_dt, C_t, dC_dt, p_t, e_t, θmin_t)

iChimera.SelfAccelerationHarmonic.Γ(vH, xH, a)
dt_dτ

DBL_DH_space = iChimera.HarmonicCoords.jBLH(xH, a)
DBL_DH = [1.0 0.0 0.0 0.0; 0.0 DBL_DH_space[1,1] DBL_DH_space[1,2] DBL_DH_space[1,3]; 0.0 DBL_DH_space[2,1] DBL_DH_space[2,2] DBL_DH_space[2,3]; 0.0 DBL_DH_space[3,1] DBL_DH_space[3,2] DBL_DH_space[3,3]]
TemporalKilling_BL = [iChimera.Kerr.KerrMetric.g_μν(rr, θθ, ϕϕ, a, i, 1) for i = 1:4]
AxialKilling_BL = [iChimera.Kerr.KerrMetric.g_μν(rr, θθ, ϕϕ, a, i, 4) for i = 1:4]

-HarmonicCoords.dot(transpose(DBL_DH) * TemporalKilling_BL, aSF_H) / dt_dτ == dE_dt
HarmonicCoords.dot(transpose(DBL_DH) * AxialKilling_BL, aSF_H) / dt_dτ == dL_dt

aSF_BL

transpose(DBL_DH) * aSF_H

aSF_BL

# update orbital constants and fluxes — function takes as argument the fluxes computed at the end of the previous geodesic (which overlaps with the start of the current geodesic piece) in order to update the fluxes using the trapezium rule
dt_flux = tF - t0
E_1, dE_dt, L_1, dL_dt, Q_1, dQ_dt, C_1, dC_dt, p_1, e_1, θmin_1 = EvolveConstants.Evolve_BL(dt_flux, a, rr, θθ, ϕϕ, dt_dτ, r_dot, θ_dot, ϕ_dot, aSF_BL, E_t, dE_dt, L_t, dL_dt, Q_t, dQ_dt, C_t, dC_dt, p_t, e_t, θmin_t)

E_t = E_1; L_t = L_1; Q_t = Q_1; C_t = C_1; p_t = p_1; e_t = e_1; θmin_t = θmin_1;
# save constants of motion
t_Fluxes[idx_save_2] = tF;
E_arr[idx_save_2] = E_t;
E_dot_arr[idx_save_2] = dE_dt;
L_arr[idx_save_2] = L_t; 
L_dot_arr[idx_save_2] = dL_dt;
C_arr[idx_save_2] = C_t;
C_dot_arr[idx_save_2] = dC_dt;
Q_arr[idx_save_2] = Q_t;
Q_dot_arr[idx_save_2] = dQ_dt;
p_arr[idx_save_2] = p_t;
e_arr[idx_save_2] = e_t;
θmin_arr[idx_save_2] = θmin_t;
idx_save_2 += 1

# save constants and fluxes
if idx_save_2 == save_every + 1
    Inspiral.save_constants!(file, save_every, t_Fluxes, E_arr, E_dot_arr, L_arr, L_dot_arr, Q_arr, Q_dot_arr, C_arr, C_dot_arr, p_arr, e_arr, θmin_arr, save_constants, save_fluxes)
    idx_save_2 = 1
    flush(file)
end

# update ODE params
zm = cos(θmin_t)^2
zp = C_t / (a^2 * (1.0-E_t^2) * zm)    # Eq. E23
ra=p_t / (1.0 - e_t); rp=p_t / (1.0 + e_t);
A = 1.0 / (1.0 - E_t^2) - (ra + rp) / 2.0    # Eq. E20
B = a^2 * C_t / ((1.0 - E_t^2) * ra * rp)    # Eq. E21
r3 = A + sqrt(A^2 - B); r4 = A - sqrt(A^2 - B);    # Eq. E19
p3 = r3 * (1.0 - e_t); p4 = r4 * (1.0 + e_t)    # Above Eq. 96
integrator.p = @SArray [a, E_t, L_t, p_t, e_t, θmin_t, p3, p4, zp, zm];
t0 = tF
# end
# end
print("Completion: 100%   \r")

# save remaining data
if idx_save_1 != 1
    @views Inspiral.save_traj!(file, idx_save_1-1, time, r, theta, phi, gamma, save_traj, save_gamma)
    @views Inspiral.save_moments!(file, idx_save_1-1, Mij2_data, Sij2_data, Mijk3_data, Sijk3_data, Mijkl4_data)
end

if idx_save_2 != 1
    @views Inspiral.save_constants!(file, idx_save_2-1, t_Fluxes, E_arr, E_dot_arr, L_arr, L_dot_arr, Q_arr, Q_dot_arr, C_arr, C_dot_arr, p_arr, e_arr, θmin_arr, save_constants, save_fluxes)
end

if JIT
    rm(sol_filename)
    println("JIT compilation run complete.")
else
    println("File created: " * sol_filename)
end
close(file)