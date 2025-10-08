module Inspiral
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
using ..SelfAccelerationCartesian
using ..EvolveConstants
using ..Waveform
using ..HarmonicCoordDerivs
using ..CoordinateDerivs
using ..MultipoleDerivs
using ..MultipoleDerivs
using ..RRPotentials

Z_1(a::Float64) = 1 + (1 - a^2 / 1.0)^(1/3) * ((1 + a)^(1/3) + (1 - a)^(1/3))
Z_2(a::Float64) = sqrt(3 * a^2 + Z_1(a)^2)
LSO_r(a::Float64) = (3 + Z_2(a) - sqrt((3 - Z_1(a)) * (3 + Z_1(a) * 2 * Z_2(a))))   # retrograde LSO
LSO_p(a::Float64) = (3 + Z_2(a) + sqrt((3 - Z_1(a)) * (3 + Z_1(a) * 2 * Z_2(a))))   # prograde LSO

"""
    compute_inspiral(args...)

Evolve inspiral with Boyer-Lindquist coordinate time parameterization and fully analyitc computation of the approximate self-force (as opposed to the Fourier fitting approach).

- `tInspiral::Float64`: total coordinate time to evolve the inspiral.
- `compute_SF::Float64`: BL time interval between self-force computations.
- `q::Float64`: mass ratio.
- `a::Float64`: black hole spin 0 < a < 1.
- `p::Float64`: initial semi-latus rectum.
- `e::Float64`: initial eccentricity.
- `θmin::Float64`: initial inclination angle.
- `sign_Lz::Int64`: sign of the z-component of the angular momentum (+1 for prograde, -1 for retrograde).
- `psi_0::Float64`: initial radial angle variable.
- `chi_0::Float64`: initial polar angle variable.
- `phi_0::Float64`: initial azimuthal angle.
- `reltol`: relative tolerance for ODE solver.
- `abstol`: absolute tolerance for ODE solver.
- `JIT::Bool`: dummy run to JIT compile function.
- `data_path::String`: path to save data.
- `lmax_mass::Int64`: maximum mass-type multipole moment l mode to include in the flux and waveform computation with 2 ≤ lmax ≤ 4
- `lmax_current::Int64` maximum current-type multipole moment l mode to include in the flux and waveform computation with 1 ≤ lmax ≤ 3 (lmax = 1 excludes any current-type moment and only up to l=3 included at this time)
- `save_traj::Bool`: whether to save the trajectory data.
- `save_constants::Bool`: whether to save the constants of motion.
- `save_fluxes::Bool`: whether to save the fluxes.
- `save_gamma::Bool`: whether to save the Lorentz factor.
- `dt_save::Float64`: time interval between saving trajectory data.
- `save_every::Int64`: number of points in each chunk of data when saving to file.
"""

function compute_inspiral(a::Float64, p::Float64, e::Float64, θmin::Float64, sign_Lz::Int64, q::Float64, psi_0::Float64, chi_0::Float64, phi_0::Float64, compute_SF::Float64, tInspiral::Float64, dt_save::Float64, save_every::Int64, reltol::Float64=1e-14, abstol::Float64=1e-14, OnePN::Float64=1.0, TwoPN::Float64=1.0, TwoPointFivePN::Float64=1.0, coordinates::String="cartesian"; data_path::String="Data/", JIT::Bool=false, lmax_mass_fluxes::Int64, lmax_current_fluxes::Int64, save_traj::Bool, save_constants::Bool, save_fluxes::Bool, save_gamma::Bool, maxiters::Int64=Int(1e8))
    # create solution file
    sol_filename=solution_fname(a, p, e, θmin, q, psi_0, chi_0, phi_0, lmax_mass_fluxes, lmax_current_fluxes, coordinates, data_path)
    
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

    Mij5 = @MArray zeros(3, 3)
    Mij6 = @MArray zeros(3, 3)
    Mij7 = @MArray zeros(3, 3)
    Mij8 = @MArray zeros(3, 3)
    dxmMij5 = @MArray zeros(3, 3, 3)
    dxmMij6 = @MArray zeros(3, 3, 3)
    dxmMij7 = @MArray zeros(3, 3, 3)

    Mijk7 = @MArray zeros(3, 3, 3)
    Mijk8 = @MArray zeros(3, 3, 3)
    dxmMijk7 = @MArray zeros(3, 3, 3, 3)

    Sij5 = @MArray zeros(3, 3)
    Sij6 = @MArray zeros(3, 3)
    dxmSij5 = @MArray zeros(3, 3, 3)
    
    # initialize derivative arrays
    xBL = @MArray zeros(3); vBL = @MArray zeros(3); aBL = @MArray zeros(3);
    
    dxBL_dt = @MArray zeros(3); d2xBL_dt = @MArray zeros(3); d3xBL_dt = @MArray zeros(3); d4xBL_dt = @MArray zeros(3);
    d5xBL_dt = @MArray zeros(3); d6xBL_dt = @MArray zeros(3); d7xBL_dt = @MArray zeros(3); d8xBL_dt = @MArray zeros(3); d9xBL_dt = @MArray zeros(3);

    dx_dλ = @MArray zeros(3); d2x_dλ = @MArray zeros(3); d3x_dλ = @MArray zeros(3); d4x_dλ = @MArray zeros(3);
    d5x_dλ = @MArray zeros(3); d6x_dλ = @MArray zeros(3); d7x_dλ = @MArray zeros(3); d8x_dλ = @MArray zeros(3); d9x_dλ = @MArray zeros(3);

    xH = @MArray zeros(3); dxH_dt = @MArray zeros(3); d2xH_dt = @MArray zeros(3); d3xH_dt = @MArray zeros(3); d4xH_dt = @MArray zeros(3);
    d5xH_dt = @MArray zeros(3); d6xH_dt = @MArray zeros(3); d7xH_dt = @MArray zeros(3); d8xH_dt = @MArray zeros(3); d9xH_dt = @MArray zeros(3);

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
    tt, rr, θθ, ϕϕ, r_dot, θ_dot, ϕ_dot, r_ddot, θ_ddot, ϕ_ddot, dt_dτ, psi, chi = compute_geodesic_arrays(integrator, a, E_t, L_t, p_t, e_t, θmin_t, p3, p4, zp, zm)

    # COMPUTE BL COORDINATE DERIVATIVES
    xBL[1] = rr; xBL[2] = θθ; xBL[3] = ϕϕ;
    vBL[1] = r_dot; vBL[2] = θ_dot; vBL[3] = ϕ_dot;
    aBL[1] = r_ddot; aBL[2] = θ_ddot; aBL[3] = ϕ_ddot;
    CoordinateDerivs.ComputeDerivs!(xBL, sign(vBL[1]), sign(vBL[2]), dxBL_dt, d2xBL_dt, d3xBL_dt, d4xBL_dt, d5xBL_dt, d6xBL_dt, d7xBL_dt, d8xBL_dt, d9xBL_dt, dx_dλ, d2x_dλ, d3x_dλ, d4x_dλ, d5x_dλ, d6x_dλ, d7x_dλ, d8x_dλ, d9x_dλ, a, E_t, L_t, C_t);

    # COMPUTE HARMONIC COORDINATE DERIVATIVES
    HarmonicCoordDerivs.compute_harmonic_derivs!(xBL, dxBL_dt, d2xBL_dt, d3xBL_dt, d4xBL_dt, d5xBL_dt, d6xBL_dt, d7xBL_dt, d8xBL_dt, d9xBL_dt, xH, dxH_dt, d2xH_dt, d3xH_dt, d4xH_dt, d5xH_dt, d6xH_dt, d7xH_dt, d8xH_dt, d9xH_dt, a);

    # UPDATE MULTIPOLE MOMENTS AND TRAJECTORY ARRAYS
    MultipoleDerivs.compute_WF_moments!(q, Mij2, Mijk3, Mijkl4, Sij2, Sijk3, xH, dxH_dt, d2xH_dt, d3xH_dt, d4xH_dt); # note that all implemented moments are calculated for the waveform. In the post-processing step, one can choose which moments to include
    update_waveform_arrays!(idx_save_1, Mij2_data, Sij2_data, Mijk3_data, Sijk3_data, Mijkl4_data, Mij2, Sij2, Mijk3, Sijk3, Mijkl4)
    update_trajectory_arrays!(integrator, idx_save_1, time, r, theta, phi, gamma, a, E_t, L_t, p_t, e_t, θmin_t, p3, p4, zp, zm, save_traj, save_gamma)
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
    SF_step = false

    t0 = 0.0
    while integrator.t < tInspiral
        print("Completion: $(round(100 * t0/tInspiral; digits=5))%   \r")
        flush(stdout)
        
        # successively step the integrator to the next time in either times_SF or times_WF, whichever is smaller
        if length(times_SF) == 0 && length(times_WF) == 0
            break
        elseif length(times_WF) == 0
            tF = times_SF[1]
            popfirst!(times_SF)
            SF_step = true
            WF_step = false
        elseif length(times_SF) == 0
            tF = times_WF[1]
            popfirst!(times_WF)
            WF_step = true
            SF_step = false
        elseif times_SF[1] < times_WF[1]
            tF = times_SF[1]
            popfirst!(times_SF)
            SF_step = true
            WF_step = false
        else
            tF = times_WF[1]
            popfirst!(times_WF)
            WF_step = true
            SF_step = false
        end

        time_step = tF - integrator.t
        step!(integrator, time_step, true)


        ## COMPUTE HIGH-ORDER TIME DERIVATIVES OF THE TRAJECTORY ##
        # EXTRACT GEODESIC TRAJECTORY
        tt, rr, θθ, ϕϕ, r_dot, θ_dot, ϕ_dot, r_ddot, θ_ddot, ϕ_ddot, dt_dτ, psi, chi = compute_geodesic_arrays(integrator, a, E_t, L_t, p_t, e_t, θmin_t, p3, p4, zp, zm)

        # COMPUTE BL COORDINATE DERIVATIVES
        xBL[1] = rr; xBL[2] = θθ; xBL[3] = ϕϕ;
        vBL[1] = r_dot; vBL[2] = θ_dot; vBL[3] = ϕ_dot;
        aBL[1] = r_ddot; aBL[2] = θ_ddot; aBL[3] = ϕ_ddot;
        CoordinateDerivs.ComputeDerivs!(xBL, sign(vBL[1]), sign(vBL[2]), dxBL_dt, d2xBL_dt, d3xBL_dt, d4xBL_dt, d5xBL_dt, d6xBL_dt, d7xBL_dt, d8xBL_dt, d9xBL_dt, dx_dλ, d2x_dλ, d3x_dλ, d4x_dλ, d5x_dλ, d6x_dλ, d7x_dλ, d8x_dλ, d9x_dλ, a, E_t, L_t, C_t);

        # COMPUTE HARMONIC COORDINATE DERIVATIVES
        HarmonicCoordDerivs.compute_harmonic_derivs!(xBL, dxBL_dt, d2xBL_dt, d3xBL_dt, d4xBL_dt, d5xBL_dt, d6xBL_dt, d7xBL_dt, d8xBL_dt, d9xBL_dt, xH, dxH_dt, d2xH_dt, d3xH_dt, d4xH_dt, d5xH_dt, d6xH_dt, d7xH_dt, d8xH_dt, d9xH_dt, a);

        ## COMPUTE THE FLUXES AND UPDATE ORBITAL CONSTANTS OR COMPUTE WAVEFORM MOMENTS, DEPENDING ON WHETHER IT IS A SELF-FORCE OR WAVEFORM STEP ##
        if WF_step
            if idx_save_1 == save_every + 1
                Inspiral.save_traj!(file, save_every, time, r, theta, phi, gamma, save_traj, save_gamma)
                Inspiral.save_moments!(file, save_every, Mij2_data, Sij2_data, Mijk3_data, Sijk3_data, Mijkl4_data)
                idx_save_1 = 1
                flush(file)
            end

            MultipoleDerivs.compute_WF_moments!(q, Mij2, Mijk3, Mijkl4, Sij2, Sijk3, xH, dxH_dt, d2xH_dt, d3xH_dt, d4xH_dt); # note that all implemented moments are calculated for the waveform. In the post-processing step, one can choose which moments to include
            update_waveform_arrays!(idx_save_1, Mij2_data, Sij2_data, Mijk3_data, Sijk3_data, Mijkl4_data, Mij2, Sij2, Mijk3, Sijk3, Mijkl4)
            update_trajectory_arrays!(integrator, idx_save_1, time, r, theta, phi, gamma, a, E_t, L_t, p_t, e_t, θmin_t, p3, p4, zp, zm, save_traj, save_gamma)
            idx_save_1 += 1
        else            
            ##### COMPUTE SELF-FORCE AND SELECT OUT DESIRED MULTIPOLE MOMENTS TO CONTRIBUTE TO SELF-FORCE COMPUTATION (MASS QUADRUPOLE ALWAYS INCLUDED. lmax_mass_fluxes=3 INCLUDES MASS OCTUPOLE, lmax_current_fluxes=2 INCLUDES CURRENT QUADRUPOLE) #####
            MultipoleDerivs.compute_SF_moments!(q, Mij5, Mij6, Mij7, Mij8, dxmMij5, dxmMij6, dxmMij7, Mijk7, Mijk8, dxmMijk7, Sij5, Sij6, dxmSij5, xH, dxH_dt, d2xH_dt, d3xH_dt, d4xH_dt, d5xH_dt, d6xH_dt, d7xH_dt, d8xH_dt, d9xH_dt, OnePN, TwoPN, TwoPointFivePN)

            if lmax_mass_fluxes != 3
                Mijk7 .= 0.0;
                Mijk8 .= 0.0;
                dxmMijk7 .= 0.0;
            end

            if lmax_current_fluxes != 2
                Sij5 .= 0.0;
                Sij6 .= 0.0;
                dxmSij5 .= 0.0;
            end

            Vrr, ∂Vrr_∂t = RRPotentials.compute_RR_potentials!(Virr, ∂Vrr_∂a, ∂Virr_∂t, ∂Virr_∂a, xH, dxH_dt, Mij5, Mij6, Mij7, Mij8, dxmMij5, dxmMij6, dxmMij7, Mijk7, Mijk8, dxmMijk7, Sij5, Sij6, dxmSij5);

            vH .= dxH_dt
            aH .= d2xH_dt
            rH = SelfAccelerationHarmonic.norm_3d(xH);
            v = SelfAccelerationHarmonic.norm_3d(vH);
            if coordinates == "harmonic"
                SelfAccelerationHarmonic.aRRα(aSF_H, aSF_BL, xH, v, vH, xBL, rH, a, Vrr, ∂Vrr_∂t, Virr, ∂Vrr_∂a, ∂Virr_∂t, ∂Virr_∂a)
            elseif coordinates == "cartesian"
                SelfAccelerationCartesian.aRRα(aSF_H, aSF_BL, xH, v, vH, xBL, rH, a, Vrr, ∂Vrr_∂t, Virr, ∂Vrr_∂a, ∂Virr_∂t, ∂Virr_∂a)
            end

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
        end
    end
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
end

function compute_geodesic_arrays(integrator, a::Float64, E::Float64, L::Float64, p::Float64, e::Float64, θmin::Float64, p3::Float64, p4::Float64, zp::Float64, zm::Float64)
    # deconstruct solution
    t = integrator.t;
    psi = integrator.u[1];
    chi = mod.(integrator.u[2], 2π);
    ϕ = integrator.u[3];

    # compute time derivatives
    psi_dot = BLTimeGeodesics.psi_dot(psi, chi, ϕ, a, E, L, p, e, θmin, p3, p4, zp, zm)
    chi_dot = BLTimeGeodesics.chi_dot(psi, chi, ϕ, a, E, L, p, e, θmin, p3, p4, zp, zm)
    ϕ_dot = BLTimeGeodesics.phi_dot(psi, chi, ϕ, a, E, L, p, e, θmin, p3, p4, zp, zm)

    # compute BL coordinates t, r, θ and their time derivatives
    r = BLTimeGeodesics.r(psi, p, e)
    θ = acos((π/2<chi<1.5π) ? -sqrt(BLTimeGeodesics.z(chi, θmin)) : sqrt(BLTimeGeodesics.z(chi, θmin)))
    r_dot = BLTimeGeodesics.dr_dt(psi_dot, psi, p, e);
    θ_dot = BLTimeGeodesics.dθ_dt(chi_dot, chi, θ, θmin);
    v = [r_dot, θ_dot, ϕ_dot];
    dt_dτ = BLTimeGeodesics.Γ(r, θ, ϕ, v, a)

    # substitute solution back into geodesic equation to find second derivatives of BL coordinates (wrt t)
    r_ddot = BLTimeGeodesics.dr2_dt2(r, θ, ϕ, r_dot, θ_dot, ϕ_dot, a)
    θ_ddot = BLTimeGeodesics.dθ2_dt2(r, θ, ϕ, r_dot, θ_dot, ϕ_dot, a)
    ϕ_ddot = BLTimeGeodesics.dϕ2_dt2(r, θ, ϕ, r_dot, θ_dot, ϕ_dot, a)

    return t, r, θ, ϕ, r_dot, θ_dot, ϕ_dot, r_ddot, θ_ddot, ϕ_ddot, dt_dτ, psi, chi
end

function update_waveform_arrays!(idx_save::Int64, Mij2_data::AbstractArray, Sij2_data::AbstractArray, Mijk3_data::AbstractArray, Sijk3_data::AbstractArray, Mijkl4_data::AbstractArray, Mij2_data_wf_temp::AbstractArray, Sij2_data_wf_temp::AbstractArray, Mijk3_data_wf_temp::AbstractArray, Sijk3_data_wf_temp::AbstractArray, Mijkl4_data_wf_temp::AbstractArray)
    @inbounds for i = 1:3, j = 1:3
        Mij2_data[i, j][idx_save] = Mij2_data_wf_temp[i, j];
        Sij2_data[i, j][idx_save] = Sij2_data_wf_temp[i, j];
        @inbounds for k = 1:3
            Mijk3_data[i, j, k][idx_save] = Mijk3_data_wf_temp[i, j, k];
            Sijk3_data[i, j, k][idx_save] = Sijk3_data_wf_temp[i, j, k];
            @inbounds for l = 1:3
                Mijkl4_data[i, j, k, l][idx_save] = Mijkl4_data_wf_temp[i, j, k, l];
            end
        end
    end
end

function update_trajectory_arrays!(integrator, idx_save::Int64, time::Vector{Float64}, r::Vector{Float64}, theta::Vector{Float64}, phi::Vector{Float64}, gamma::Vector{Float64}, a::Float64, E::Float64, L::Float64, p::Float64, e::Float64, θmin::Float64, p3::Float64, p4::Float64, zp::Float64, zm::Float64, save_traj::Bool, save_gamma::Bool)
    t = integrator.t;
    psi = integrator.u[1];
    chi = mod.(integrator.u[2], 2π);
    ϕ = integrator.u[3];

    time[idx_save] = t

    if save_traj || save_gamma
        tt, rr, θθ, ϕϕ, r_dot, θ_dot, ϕ_dot, r_ddot, θ_ddot, ϕ_ddot, dt_dτ, psi, chi = compute_geodesic_arrays(integrator, a, E, L, p, e, θmin, p3, p4, zp, zm)
    end

    if save_traj
        r[idx_save] = rr
        theta[idx_save] = θθ
        phi[idx_save] = ϕϕ
    end

    if save_gamma
        gamma[idx_save] = dt_dτ
    end
end


# evolve inspiral along one piecewise geodesic
function evolve_inspiral!(integrator, h::Number, tt::Vector{<:Number}, dt_dτ::Vector{<:Number}, rr::Vector{<:Number}, r_dot::Vector{<:Number}, r_ddot::Vector{<:Number}, θθ::Vector{<:Number}, θ_dot::Vector{<:Number}, θ_ddot::Vector{<:Number}, 
    ϕϕ::Vector{<:Number}, ϕ_dot::Vector{<:Number}, ϕ_ddot::Vector{<:Number}, dt_dλλ::Vector{<:Number})
    a, E, L, p, e, θmin, p3, p4, zp, zm = integrator.p
    track_num_steps = 0
    @inbounds for i = 1:length(tt)
        track_num_steps += 1
        compute_BL_coords_traj!(integrator, i, λλ, tt, dt_dτ, rr, r_dot, r_ddot, θθ, θ_dot, θ_ddot, ϕϕ, ϕ_dot, ϕ_ddot, dt_dλλ, a, E, L, p, e, θmin, p3, p4, zp, zm)
        step!(integrator, h, true)
    end
    if track_num_steps != length(λλ)
        throw(ArgumentError("Length of λλ array does not match the number ($(i)) of steps taken"))
    end
end

function solution_fname(a::Float64, p::Float64, e::Float64, θmin::Float64, q::Float64, psi_0::Float64, chi_0::Float64, phi_0::Float64, lmax_mass_fluxes::Int64, lmax_current_fluxes::Int64, coordinates::String, data_path::String)
    return data_path * "EMRI_sol_a_$(a)_p_$(p)_e_$(e)_θmin_$(round(θmin, sigdigits=3))_q_$(q)_psi0_$(round(psi_0, sigdigits=3))_chi0_$(round(chi_0, sigdigits=3))_phi0_$(round(phi_0, sigdigits=3))_BL_time_lmax_mass_fluxes_$(lmax_mass_fluxes)_lmax_current_fluxes_$(lmax_current_fluxes)_" * coordinates * ".h5"
end

function waveform_fname(a::Float64, p::Float64, e::Float64, θmin::Float64, q::Float64, psi_0::Float64, chi_0::Float64, phi_0::Float64, obs_distance::Float64, ThetaSource::Float64, PhiSource::Float64, ThetaKerr::Float64, PhiKerr::Float64, lmax_mass_fluxes::Int64, lmax_current_fluxes::Int64, lmax_mass_waveform::Int64, lmax_current_waveform::Int64, coordinates::String, data_path::String)
    return data_path * "Waveform_a_$(a)_p_$(p)_e_$(e)_θmin_$(round(θmin, sigdigits=3))_q_$(q)_psi0_$(round(psi_0, sigdigits=3))_chi0_$(round(chi_0, sigdigits=3))_phi0_$(round(phi_0, sigdigits=3))_obsDist_$(round(obs_distance, sigdigits=3))_ThetaS_$(round(ThetaSource, sigdigits=3))_PhiS_$(round(PhiSource, sigdigits=3))_ThetaK_$(round(ThetaKerr, sigdigits=3))_PhiK_$(round(PhiKerr, sigdigits=3))_BL_time_lmax_mass_fluxes_$(lmax_mass_fluxes)_lmax_current_fluxes_$(lmax_current_fluxes)_lmax_mass_waveform_$(lmax_mass_waveform)_lmax_current_waveform_$(lmax_current_waveform)_" * coordinates * ".h5"
end

function waveform_fname(a::Float64, p::Float64, e::Float64, θmin::Float64, q::Float64, psi_0::Float64, chi_0::Float64, phi_0::Float64, obs_distance::Float64, ThetaObs::Float64, PhiObs::Float64, lmax_mass_fluxes::Int64, lmax_current_fluxes::Int64, lmax_mass_waveform::Int64, lmax_current_waveform::Int64, coordinates::String, data_path::String)
    return data_path * "Waveform_a_$(a)_p_$(p)_e_$(e)_θmin_$(round(θmin, sigdigits=3))_q_$(q)_psi0_$(round(psi_0, sigdigits=3))_chi0_$(round(chi_0, sigdigits=3))_phi0_$(round(phi_0, sigdigits=3))_obsDist_$(round(obs_distance, sigdigits=3))_ThetaObs_$(round(ThetaObs, sigdigits=3))_PhiObs_$(round(PhiObs, sigdigits=3))_BL_time_lmax_mass_fluxes_$(lmax_mass_fluxes)_lmax_current_fluxes_$(lmax_current_fluxes)_lmax_mass_waveform_$(lmax_mass_waveform)_lmax_current_waveform_$(lmax_current_waveform)_" * coordinates * ".h5"
end

function load_trajectory(a::Float64, p::Float64, e::Float64, θmin::Float64, q::Float64, psi_0::Float64, chi_0::Float64, phi_0::Float64, lmax_mass_fluxes::Int64, lmax_current_fluxes::Int64, coordinates::String, data_path::String)
    sol_filename=solution_fname(a, p, e, θmin, q, psi_0, chi_0, phi_0, lmax_mass_fluxes, lmax_current_fluxes, coordinates, data_path)
    return Inspiral.load_trajectory(sol_filename)
end

function load_constants_of_motion(a::Float64, p::Float64, e::Float64, θmin::Float64, q::Float64, psi_0::Float64, chi_0::Float64, phi_0::Float64, lmax_mass_fluxes::Int64, lmax_current_fluxes::Int64, coordinates::String, data_path::String)
    sol_filename=solution_fname(a, p, e, θmin, q, psi_0, chi_0, phi_0, lmax_mass_fluxes, lmax_current_fluxes, coordinates, data_path)
    return Inspiral.load_constants_of_motion(sol_filename)
end

function compute_waveform(obs_distance::Float64, ThetaSource::Float64, PhiSource::Float64, ThetaKerr::Float64, PhiKerr::Float64, a::Float64, p::Float64, e::Float64, θmin::Float64, q::Float64, psi_0::Float64, chi_0::Float64, phi_0::Float64, lmax_mass_fluxes::Int64, lmax_current_fluxes::Int64, lmax_mass_waveform::Int64, lmax_current_waveform::Int64, coordinates::String,  data_path::String)
    # load waveform multipole moments
    sol_filename=solution_fname(a, p, e, θmin, q, psi_0, chi_0, phi_0, lmax_mass_fluxes, lmax_current_fluxes, coordinates, data_path)
    t, Mij2, Mijk3, Mijkl4, Sij2, Sijk3 = Inspiral.load_waveform_moments(sol_filename, lmax_mass_waveform, lmax_current_waveform)
    num_points = length(Mij2[1, 1]);
    h_plus, h_cross = Waveform.compute_wave_polarizations(num_points, obs_distance, deg2rad(ThetaSource), deg2rad(PhiSource), deg2rad(ThetaKerr), deg2rad(PhiKerr), Mij2, Mijk3, Mijkl4, Sij2, Sijk3, q)

    # save waveform to file
    wave_filename=waveform_fname(a, p, e, θmin, q, psi_0, chi_0, phi_0, obs_distance, ThetaSource, PhiSource, ThetaKerr, PhiKerr, lmax_mass_fluxes, lmax_current_fluxes, lmax_mass_waveform, lmax_current_waveform, coordinates, data_path)
    h5open(wave_filename, "w") do file
        file["t"] = t
        file["hplus"] = h_plus
        file["hcross"] = h_cross
    end
    println("File created: " * wave_filename)
end

function compute_waveform(obs_distance::Float64, ThetaObs::Float64, PhiObs::Float64, a::Float64, p::Float64, e::Float64, θmin::Float64, q::Float64, psi_0::Float64, chi_0::Float64, phi_0::Float64, lmax_mass_fluxes::Int64, lmax_current_fluxes::Int64, lmax_mass_waveform::Int64, lmax_current_waveform::Int64, coordinates::String,  data_path::String)
    # load waveform multipole moments
    sol_filename=solution_fname(a, p, e, θmin, q, psi_0, chi_0, phi_0, lmax_mass_fluxes, lmax_current_fluxes, coordinates, data_path)
    t, Mij2, Mijk3, Mijkl4, Sij2, Sijk3 = Inspiral.load_waveform_moments(sol_filename, lmax_mass_waveform, lmax_current_waveform)
    num_points = length(Mij2[1, 1]);
    h_plus, h_cross = Waveform.compute_wave_polarizations(num_points, obs_distance, deg2rad(ThetaObs), deg2rad(PhiObs), Mij2, Mijk3, Mijkl4, Sij2, Sijk3, q)

    # save waveform to file
    wave_filename=waveform_fname(a, p, e, θmin, q, psi_0, chi_0, phi_0, obs_distance, ThetaObs, PhiObs, lmax_mass_fluxes, lmax_current_fluxes, lmax_mass_waveform, lmax_current_waveform, coordinates, data_path)
    h5open(wave_filename, "w") do file
        file["t"] = t
        file["hplus"] = h_plus
        file["hcross"] = h_cross
    end
    println("File created: " * wave_filename)
end

function load_waveform(obs_distance::Float64, ThetaSource::Float64, PhiSource::Float64, ThetaKerr::Float64, PhiKerr::Float64, a::Float64, p::Float64, e::Float64, θmin::Float64, q::Float64, psi_0::Float64, chi_0::Float64, phi_0::Float64, lmax_mass_fluxes::Int64, lmax_current_fluxes::Int64, lmax_mass_waveform::Int64, lmax_current_waveform::Int64, coordinates::String,  data_path::String)
    # save waveform to file
    wave_filename=waveform_fname(a, p, e, θmin, q, psi_0, chi_0, phi_0, obs_distance, ThetaSource, PhiSource, ThetaKerr, PhiKerr, lmax_mass_fluxes, lmax_current_fluxes, lmax_mass_waveform, lmax_current_waveform, coordinates, data_path)
    file = h5open(wave_filename, "r")
    t = file["t"][:]
    h_plus = file["hplus"][:]
    h_cross = file["hcross"][:]
    close(file)
    return t, h_plus, h_cross    
end

function load_waveform(obs_distance::Float64, ThetaObs::Float64, PhiObs::Float64, a::Float64, p::Float64, e::Float64, θmin::Float64, q::Float64, psi_0::Float64, chi_0::Float64, phi_0::Float64, lmax_mass_fluxes::Int64, lmax_current_fluxes::Int64, lmax_mass_waveform::Int64, lmax_current_waveform::Int64, coordinates::String,  data_path::String)
    # save waveform to file
    wave_filename=waveform_fname(a, p, e, θmin, q, psi_0, chi_0, phi_0, obs_distance, ThetaObs, PhiObs, lmax_mass_fluxes, lmax_current_fluxes, lmax_mass_waveform, lmax_current_waveform, coordinates, data_path)
    file = h5open(wave_filename, "r")
    t = file["t"][:]
    h_plus = file["hplus"][:]
    h_cross = file["hcross"][:]
    close(file)
    return t, h_plus, h_cross    
end

# useful for dummy runs (e.g., for resonances to estimate the duration of time needed by computing the time derivative of the fundamental frequencies)
function delete_EMRI_data(a::Float64, p::Float64, e::Float64, θmin::Float64, q::Float64, psi_0::Float64, chi_0::Float64, phi_0::Float64, lmax_mass_fluxes::Int64, lmax_current_fluxes::Int64, coordinates::String, data_path::String)
    sol_filename=solution_fname(a, p, e, θmin, q, psi_0, chi_0, phi_0, lmax_mass_fluxes, lmax_current_fluxes, coordinates, data_path)
    rm(sol_filename)
end

function initialize_solution_file!(file::HDF5.File, chunk_size::Int64, save_traj::Bool, save_constants::Bool, save_fluxes::Bool, save_gamma::Bool)
    create_dataset!(file, "", "t", Float64, chunk_size);

    traj_group_name = "Trajectory"
    create_file_group!(file, traj_group_name);
    
    if save_traj
        create_dataset!(file, traj_group_name, "r", Float64, chunk_size);
        create_dataset!(file, traj_group_name, "theta", Float64, chunk_size);
        create_dataset!(file, traj_group_name, "phi", Float64, chunk_size);
    end

    if save_gamma
        create_dataset!(file, traj_group_name, "Gamma", Float64, chunk_size);
    end

    constants_group_name = "ConstantsOfMotion"
    create_file_group!(file, constants_group_name);
    if save_constants || save_fluxes
        create_dataset!(file, constants_group_name, "t", Float64, chunk_size);
    end

    if save_constants
        create_dataset!(file, constants_group_name, "Energy", Float64, chunk_size);
        create_dataset!(file, constants_group_name, "AngularMomentum", Float64, chunk_size);
        create_dataset!(file, constants_group_name, "CarterConstant", Float64, chunk_size);
        create_dataset!(file, constants_group_name, "AltCarterConstant", Float64, chunk_size);
        create_dataset!(file, constants_group_name, "p", Float64, chunk_size);
        create_dataset!(file, constants_group_name, "eccentricity", Float64, chunk_size);
        create_dataset!(file, constants_group_name, "theta_min", Float64, chunk_size);
    end

    if save_fluxes
        create_dataset!(file, constants_group_name, "Edot", Float64, chunk_size);
        create_dataset!(file, constants_group_name, "Ldot", Float64, chunk_size);
        create_dataset!(file, constants_group_name, "Qdot", Float64, chunk_size);
        create_dataset!(file, constants_group_name, "Cdot", Float64, chunk_size);
    end

    ## INDEPENDENT WAVEFORM MOMENT COMPONENTS ##
    wave_group_name = "WaveformMoments"
    create_file_group!(file, wave_group_name);

    # mass and current quadrupole second time derivs
    create_dataset!(file, wave_group_name, "Mij11_2", Float64, chunk_size); create_dataset!(file, wave_group_name, "Mij12_2", Float64, chunk_size);  create_dataset!(file, wave_group_name, "Mij13_2", Float64, chunk_size);
    create_dataset!(file, wave_group_name, "Mij22_2", Float64, chunk_size); create_dataset!(file, wave_group_name, "Mij23_2", Float64, chunk_size);  create_dataset!(file, wave_group_name, "Mij33_2", Float64, chunk_size);

    create_dataset!(file, wave_group_name, "Sij11_2", Float64, chunk_size); create_dataset!(file, wave_group_name, "Sij12_2", Float64, chunk_size);  create_dataset!(file, wave_group_name, "Sij13_2", Float64, chunk_size);
    create_dataset!(file, wave_group_name, "Sij22_2", Float64, chunk_size); create_dataset!(file, wave_group_name, "Sij23_2", Float64, chunk_size);  create_dataset!(file, wave_group_name, "Sij33_2", Float64, chunk_size);

    # mass and current octupole third time derivs
    create_dataset!(file, wave_group_name, "Mijk111_3", Float64, chunk_size); create_dataset!(file, wave_group_name, "Mijk112_3", Float64, chunk_size); create_dataset!(file, wave_group_name, "Mijk122_3", Float64, chunk_size); 
    create_dataset!(file, wave_group_name, "Mijk113_3", Float64, chunk_size); create_dataset!(file, wave_group_name, "Mijk133_3", Float64, chunk_size); create_dataset!(file, wave_group_name, "Mijk123_3", Float64, chunk_size); 
    create_dataset!(file, wave_group_name, "Mijk222_3", Float64, chunk_size); create_dataset!(file, wave_group_name, "Mijk223_3", Float64, chunk_size); create_dataset!(file, wave_group_name, "Mijk233_3", Float64, chunk_size); 
    create_dataset!(file, wave_group_name, "Mijk333_3", Float64, chunk_size);

    create_dataset!(file, wave_group_name, "Sijk111_3", Float64, chunk_size); create_dataset!(file, wave_group_name, "Sijk112_3", Float64, chunk_size); create_dataset!(file, wave_group_name, "Sijk122_3", Float64, chunk_size); 
    create_dataset!(file, wave_group_name, "Sijk113_3", Float64, chunk_size); create_dataset!(file, wave_group_name, "Sijk133_3", Float64, chunk_size); create_dataset!(file, wave_group_name, "Sijk123_3", Float64, chunk_size); 
    create_dataset!(file, wave_group_name, "Sijk222_3", Float64, chunk_size); create_dataset!(file, wave_group_name, "Sijk223_3", Float64, chunk_size); create_dataset!(file, wave_group_name, "Sijk233_3", Float64, chunk_size); 
    create_dataset!(file, wave_group_name, "Sijk333_3", Float64, chunk_size);


    # mass hexadecapole fourth time deriv
    create_dataset!(file, wave_group_name, "Mijkl1111_4", Float64, chunk_size); create_dataset!(file, wave_group_name, "Mijkl1112_4", Float64, chunk_size); create_dataset!(file, wave_group_name, "Mijkl1122_4", Float64, chunk_size);
    create_dataset!(file, wave_group_name, "Mijkl1222_4", Float64, chunk_size); create_dataset!(file, wave_group_name, "Mijkl1113_4", Float64, chunk_size); create_dataset!(file, wave_group_name, "Mijkl1133_4", Float64, chunk_size);
    create_dataset!(file, wave_group_name, "Mijkl1333_4", Float64, chunk_size); create_dataset!(file, wave_group_name, "Mijkl1123_4", Float64, chunk_size); create_dataset!(file, wave_group_name, "Mijkl1223_4", Float64, chunk_size);
    create_dataset!(file, wave_group_name, "Mijkl1233_4", Float64, chunk_size); create_dataset!(file, wave_group_name, "Mijkl2222_4", Float64, chunk_size); create_dataset!(file, wave_group_name, "Mijkl2223_4", Float64, chunk_size);
    create_dataset!(file, wave_group_name, "Mijkl2233_4", Float64, chunk_size); create_dataset!(file, wave_group_name, "Mijkl2333_4", Float64, chunk_size); create_dataset!(file, wave_group_name, "Mijkl3333_4", Float64, chunk_size);

    return file
end

@views function save_traj!(file::HDF5.File, chunk_size::Int64, t::Vector{Float64}, r::Vector{Float64}, θ::Vector{Float64}, ϕ::Vector{Float64}, dt_dτ::Vector{Float64}, save_traj::Bool, save_gamma::Bool)
    append_data!(file, "", "t", t[1:chunk_size], chunk_size);
    traj_group_name = "Trajectory"

    if save_traj
        append_data!(file, traj_group_name, "r", r[1:chunk_size], chunk_size);
        append_data!(file, traj_group_name, "theta", θ[1:chunk_size], chunk_size);
        append_data!(file, traj_group_name, "phi", ϕ[1:chunk_size], chunk_size);
    end

    if save_gamma
        append_data!(file, traj_group_name, "Gamma", dt_dτ[1:chunk_size], chunk_size);
    end
end

@views function save_moments!(file::HDF5.File, chunk_size::Int64, Mij2::AbstractArray, Sij2::AbstractArray, Mijk3::AbstractArray, Sijk3::AbstractArray, Mijkl4::AbstractArray)

    ## INDEPENDENT WAVEFORM MOMENT COMPONENTS ##
    wave_group_name = "WaveformMoments"
    # mass and current quadrupole second time derivs
    append_data!(file, wave_group_name, "Mij11_2", Mij2[1,1][1:chunk_size], chunk_size); append_data!(file, wave_group_name, "Mij12_2", Mij2[1,2][1:chunk_size], chunk_size);  append_data!(file, wave_group_name, "Mij13_2", Mij2[1,3][1:chunk_size], chunk_size);
    append_data!(file, wave_group_name, "Mij22_2", Mij2[2,2][1:chunk_size], chunk_size); append_data!(file, wave_group_name, "Mij23_2", Mij2[2,3][1:chunk_size], chunk_size);  append_data!(file, wave_group_name, "Mij33_2", Mij2[3,3][1:chunk_size], chunk_size);

    append_data!(file, wave_group_name, "Sij11_2", Sij2[1,1][1:chunk_size], chunk_size); append_data!(file, wave_group_name, "Sij12_2", Sij2[1,2][1:chunk_size], chunk_size);  append_data!(file, wave_group_name, "Sij13_2", Sij2[1,3][1:chunk_size], chunk_size);
    append_data!(file, wave_group_name, "Sij22_2", Sij2[2,2][1:chunk_size], chunk_size); append_data!(file, wave_group_name, "Sij23_2", Sij2[2,3][1:chunk_size], chunk_size);  append_data!(file, wave_group_name, "Sij33_2", Sij2[3,3][1:chunk_size], chunk_size);


    # mass and current octupole third time derivs
    append_data!(file, wave_group_name, "Mijk111_3", Mijk3[1,1,1][1:chunk_size], chunk_size); append_data!(file, wave_group_name, "Mijk112_3", Mijk3[1,1,2][1:chunk_size], chunk_size); append_data!(file, wave_group_name, "Mijk122_3", Mijk3[1,2,2][1:chunk_size], chunk_size); 
    append_data!(file, wave_group_name, "Mijk113_3", Mijk3[1,1,3][1:chunk_size], chunk_size); append_data!(file, wave_group_name, "Mijk133_3", Mijk3[1,3,3][1:chunk_size], chunk_size); append_data!(file, wave_group_name, "Mijk123_3", Mijk3[1,2,3][1:chunk_size], chunk_size); 
    append_data!(file, wave_group_name, "Mijk222_3", Mijk3[2,2,2][1:chunk_size], chunk_size); append_data!(file, wave_group_name, "Mijk223_3", Mijk3[2,2,3][1:chunk_size], chunk_size); append_data!(file, wave_group_name, "Mijk233_3", Mijk3[2,3,3][1:chunk_size], chunk_size); 
    append_data!(file, wave_group_name, "Mijk333_3", Mijk3[3,3,3][1:chunk_size], chunk_size);



    append_data!(file, wave_group_name, "Sijk111_3", Sijk3[1,1,1][1:chunk_size], chunk_size); append_data!(file, wave_group_name, "Sijk112_3", Sijk3[1,1,2][1:chunk_size], chunk_size); append_data!(file, wave_group_name, "Sijk122_3", Sijk3[1,2,2][1:chunk_size], chunk_size); 
    append_data!(file, wave_group_name, "Sijk113_3", Sijk3[1,1,3][1:chunk_size], chunk_size); append_data!(file, wave_group_name, "Sijk133_3", Sijk3[1,3,3][1:chunk_size], chunk_size); append_data!(file, wave_group_name, "Sijk123_3", Sijk3[1,2,3][1:chunk_size], chunk_size); 
    append_data!(file, wave_group_name, "Sijk222_3", Sijk3[2,2,2][1:chunk_size], chunk_size); append_data!(file, wave_group_name, "Sijk223_3", Sijk3[2,2,3][1:chunk_size], chunk_size); append_data!(file, wave_group_name, "Sijk233_3", Sijk3[2,3,3][1:chunk_size], chunk_size); 
    append_data!(file, wave_group_name, "Sijk333_3", Sijk3[3,3,3][1:chunk_size], chunk_size);

    # mass hexadecapole fourth time deriv
    append_data!(file, wave_group_name, "Mijkl1111_4", Mijkl4[1,1,1,1][1:chunk_size], chunk_size); append_data!(file, wave_group_name, "Mijkl1112_4", Mijkl4[1,1,1,2][1:chunk_size], chunk_size); append_data!(file, wave_group_name, "Mijkl1122_4", Mijkl4[1,1,2,2][1:chunk_size], chunk_size);
    append_data!(file, wave_group_name, "Mijkl1222_4", Mijkl4[1,2,2,2][1:chunk_size], chunk_size); append_data!(file, wave_group_name, "Mijkl1113_4", Mijkl4[1,1,1,3][1:chunk_size], chunk_size); append_data!(file, wave_group_name, "Mijkl1133_4", Mijkl4[1,1,3,3][1:chunk_size], chunk_size);
    append_data!(file, wave_group_name, "Mijkl1333_4", Mijkl4[1,3,3,3][1:chunk_size], chunk_size); append_data!(file, wave_group_name, "Mijkl1123_4", Mijkl4[1,1,2,3][1:chunk_size], chunk_size); append_data!(file, wave_group_name, "Mijkl1223_4", Mijkl4[1,2,2,3][1:chunk_size], chunk_size);
    append_data!(file, wave_group_name, "Mijkl1233_4", Mijkl4[1,2,3,3][1:chunk_size], chunk_size); append_data!(file, wave_group_name, "Mijkl2222_4", Mijkl4[2,2,2,2][1:chunk_size], chunk_size); append_data!(file, wave_group_name, "Mijkl2223_4", Mijkl4[2,2,2,3][1:chunk_size], chunk_size);
    append_data!(file, wave_group_name, "Mijkl2233_4", Mijkl4[2,2,3,3][1:chunk_size], chunk_size); append_data!(file, wave_group_name, "Mijkl2333_4", Mijkl4[2,3,3,3][1:chunk_size], chunk_size); append_data!(file, wave_group_name, "Mijkl3333_4", Mijkl4[3,3,3,3][1:chunk_size], chunk_size);

end

@views function save_constants!(file::HDF5.File, chunk_size::Int64, t::Vector{Float64}, E::Vector{Float64}, dE_dt::Vector{Float64}, L::Vector{Float64}, dL_dt::Vector{Float64}, Q::Vector{Float64}, dQ_dt::Vector{Float64}, C::Vector{Float64}, dC_dt::Vector{Float64}, p::Vector{Float64}, e::Vector{Float64}, θmin::Vector{Float64}, save_constants::Bool, save_fluxes::Bool)
    constants_group_name = "ConstantsOfMotion"
    if save_constants || save_fluxes
        append_data!(file, constants_group_name, "t", t[1:chunk_size], chunk_size);
    end

    if save_constants
        append_data!(file, constants_group_name, "Energy", E[1:chunk_size], chunk_size);
        append_data!(file, constants_group_name, "AngularMomentum", L[1:chunk_size], chunk_size);
        append_data!(file, constants_group_name, "CarterConstant", C[1:chunk_size], chunk_size);
        append_data!(file, constants_group_name, "AltCarterConstant", Q[1:chunk_size], chunk_size);
        append_data!(file, constants_group_name, "p", p[1:chunk_size], chunk_size);
        append_data!(file, constants_group_name, "eccentricity", e[1:chunk_size], chunk_size);
        append_data!(file, constants_group_name, "theta_min", θmin[1:chunk_size], chunk_size);
    end

    if save_fluxes
        append_data!(file, constants_group_name, "Edot", dE_dt[1:chunk_size], chunk_size);
        append_data!(file, constants_group_name, "Ldot", dL_dt[1:chunk_size], chunk_size);
        append_data!(file, constants_group_name, "Qdot", dQ_dt[1:chunk_size], chunk_size);
        append_data!(file, constants_group_name, "Cdot", dC_dt[1:chunk_size], chunk_size);
    end
end

function save_self_acceleration!(file::HDF5.File, acc_BL::Vector{Float64}, acc_Harm::Vector{Float64})
    SF_group_name = "SelfForce"
    append_data!(file, SF_group_name, "self_acc_BL_t", acc_BL[1]);
    append_data!(file, SF_group_name, "self_acc_BL_r", acc_BL[2]);
    append_data!(file, SF_group_name, "self_acc_BL_θ", acc_BL[3]);
    append_data!(file, SF_group_name, "self_acc_BL_ϕ", acc_BL[4]);
    append_data!(file, SF_group_name, "self_acc_Harm_t", acc_Harm[1]);
    append_data!(file, SF_group_name, "self_acc_Harm_x", acc_Harm[2]);
    append_data!(file, SF_group_name, "self_acc_Harm_y", acc_Harm[3]);
    append_data!(file, SF_group_name, "self_acc_Harm_z", acc_Harm[4]);
end

function load_waveform_moments(sol_filename::String, lmax_mass::Int64, lmax_current::Int64)
    Mij2 = [Float64[] for i=1:3, j=1:3]
    Sij2 = [Float64[] for i=1:3, j=1:3]
    Mijk3 = [Float64[] for i=1:3, j=1:3, k=1:3]
    Sijk3 = [Float64[] for i=1:3, j=1:3, k=1:3]
    Mijkl4 = [Float64[] for i=1:3, j=1:3, k=1:3, l=1:3]
   
    h5f = h5open(sol_filename, "r")

    # mass and current quadrupole second time derivs
    Mij2[1,1] = h5f["WaveformMoments/Mij11_2"][:];
    Mij2[1,2] = h5f["WaveformMoments/Mij12_2"][:];
    Mij2[1,3] = h5f["WaveformMoments/Mij13_2"][:];
    Mij2[2,2] = h5f["WaveformMoments/Mij22_2"][:];
    Mij2[2,3] = h5f["WaveformMoments/Mij23_2"][:];
    Mij2[3,3] = h5f["WaveformMoments/Mij33_2"][:];


    if lmax_current >= 2
        Sij2[1,1] = h5f["WaveformMoments/Sij11_2"][:];
        Sij2[1,2] = h5f["WaveformMoments/Sij12_2"][:];
        Sij2[1,3] = h5f["WaveformMoments/Sij13_2"][:];
        Sij2[2,2] = h5f["WaveformMoments/Sij22_2"][:];
        Sij2[2,3] = h5f["WaveformMoments/Sij23_2"][:];
        Sij2[3,3] = h5f["WaveformMoments/Sij33_2"][:];
    else
        Sij2[1,1] = zeros(length(Mij2[1,1]));
        Sij2[1,2] = zeros(length(Mij2[1,2]));
        Sij2[1,3] = zeros(length(Mij2[1,3]));
        Sij2[2,2] = zeros(length(Mij2[2,2]));
        Sij2[2,3] = zeros(length(Mij2[2,3]));
        Sij2[3,3] = zeros(length(Mij2[3,3]));
    end


    # mass and current octupole third time derivs
    if lmax_mass >= 3
        Mijk3[1,1,1] = h5f["WaveformMoments/Mijk111_3"][:];
        Mijk3[1,1,2] = h5f["WaveformMoments/Mijk112_3"][:];
        Mijk3[1,2,2] = h5f["WaveformMoments/Mijk122_3"][:];
        Mijk3[1,1,3] = h5f["WaveformMoments/Mijk113_3"][:];
        Mijk3[1,3,3] = h5f["WaveformMoments/Mijk133_3"][:];
        Mijk3[1,2,3] = h5f["WaveformMoments/Mijk123_3"][:];
        Mijk3[2,2,2] = h5f["WaveformMoments/Mijk222_3"][:];
        Mijk3[2,2,3] = h5f["WaveformMoments/Mijk223_3"][:];
        Mijk3[2,3,3] = h5f["WaveformMoments/Mijk233_3"][:];
        Mijk3[3,3,3] = h5f["WaveformMoments/Mijk333_3"][:];
    else
        Mijk3[1,1,1] = zeros(length(Mij2[1,1]));
        Mijk3[1,1,2] = zeros(length(Mij2[1,2]));
        Mijk3[1,2,2] = zeros(length(Mij2[2,2]));
        Mijk3[1,1,3] = zeros(length(Mij2[1,3]));
        Mijk3[1,3,3] = zeros(length(Mij2[3,3]));
        Mijk3[1,2,3] = zeros(length(Mij2[2,3]));
        Mijk3[2,2,2] = zeros(length(Mij2[2,2]));
        Mijk3[2,2,3] = zeros(length(Mij2[2,3]));
        Mijk3[2,3,3] = zeros(length(Mij2[3,3]));
        Mijk3[3,3,3] = zeros(length(Mij2[3,3]));
    end

    if lmax_current == 3
        Sijk3[1,1,1] = h5f["WaveformMoments/Sijk111_3"][:];
        Sijk3[1,1,2] = h5f["WaveformMoments/Sijk112_3"][:];
        Sijk3[1,2,2] = h5f["WaveformMoments/Sijk122_3"][:];
        Sijk3[1,1,3] = h5f["WaveformMoments/Sijk113_3"][:];
        Sijk3[1,3,3] = h5f["WaveformMoments/Sijk133_3"][:];
        Sijk3[1,2,3] = h5f["WaveformMoments/Sijk123_3"][:];
        Sijk3[2,2,2] = h5f["WaveformMoments/Sijk222_3"][:];
        Sijk3[2,2,3] = h5f["WaveformMoments/Sijk223_3"][:];
        Sijk3[2,3,3] = h5f["WaveformMoments/Sijk233_3"][:];
        Sijk3[3,3,3] = h5f["WaveformMoments/Sijk333_3"][:];
    else
        Sijk3[1,1,1] = zeros(length(Mij2[1,1]));
        Sijk3[1,1,2] = zeros(length(Mij2[1,2]));
        Sijk3[1,2,2] = zeros(length(Mij2[2,2]));
        Sijk3[1,1,3] = zeros(length(Mij2[1,3]));
        Sijk3[1,3,3] = zeros(length(Mij2[3,3]));
        Sijk3[1,2,3] = zeros(length(Mij2[2,3]));
        Sijk3[2,2,2] = zeros(length(Mij2[2,2]));
        Sijk3[2,2,3] = zeros(length(Mij2[2,3]));
        Sijk3[2,3,3] = zeros(length(Mij2[3,3]));
        Sijk3[3,3,3] = zeros(length(Mij2[3,3]));
    end


    # mass hexadecapole fourth time deriv
    if lmax_mass == 4
        Mijkl4[1,1,1,1] = h5f["WaveformMoments/Mijkl1111_4"][:];
        Mijkl4[1,1,1,2] = h5f["WaveformMoments/Mijkl1112_4"][:];
        Mijkl4[1,1,2,2] = h5f["WaveformMoments/Mijkl1122_4"][:];
        Mijkl4[1,2,2,2] = h5f["WaveformMoments/Mijkl1222_4"][:];
        Mijkl4[1,1,1,3] = h5f["WaveformMoments/Mijkl1113_4"][:];
        Mijkl4[1,1,3,3] = h5f["WaveformMoments/Mijkl1133_4"][:];
        Mijkl4[1,3,3,3] = h5f["WaveformMoments/Mijkl1333_4"][:];
        Mijkl4[1,1,2,3] = h5f["WaveformMoments/Mijkl1123_4"][:];
        Mijkl4[1,2,2,3] = h5f["WaveformMoments/Mijkl1223_4"][:];
        Mijkl4[1,2,3,3] = h5f["WaveformMoments/Mijkl1233_4"][:];
        Mijkl4[2,2,2,2] = h5f["WaveformMoments/Mijkl2222_4"][:];
        Mijkl4[2,2,2,3] = h5f["WaveformMoments/Mijkl2223_4"][:];
        Mijkl4[2,2,3,3] = h5f["WaveformMoments/Mijkl2233_4"][:];
        Mijkl4[2,3,3,3] = h5f["WaveformMoments/Mijkl2333_4"][:];
        Mijkl4[3,3,3,3] = h5f["WaveformMoments/Mijkl3333_4"][:];
    else
        Mijkl4[1,1,1,1] = zeros(length(Mij2[1,1]));
        Mijkl4[1,1,1,2] = zeros(length(Mij2[1,2]));
        Mijkl4[1,1,2,2] = zeros(length(Mij2[2,2]));
        Mijkl4[1,2,2,2] = zeros(length(Mij2[2,2]));
        Mijkl4[1,1,1,3] = zeros(length(Mij2[1,3]));
        Mijkl4[1,1,3,3] = zeros(length(Mij2[3,3]));
        Mijkl4[1,3,3,3] = zeros(length(Mij2[3,3]));
        Mijkl4[1,1,2,3] = zeros(length(Mij2[2,3]));
        Mijkl4[1,2,2,3] = zeros(length(Mij2[2,3]));
        Mijkl4[1,2,3,3] = zeros(length(Mij2[3,3]));
        Mijkl4[2,2,2,2] = zeros(length(Mij2[2,2]));
        Mijkl4[2,2,2,3] = zeros(length(Mij2[2,3]));
        Mijkl4[2,2,3,3] = zeros(length(Mij2[3,3]));
        Mijkl4[2,3,3,3] = zeros(length(Mij2[3,3]));
        Mijkl4[3,3,3,3] = zeros(length(Mij2[3,3]));
    end

    t = h5f["t"][:];
    close(h5f)

    # symmetrize 
    SymmetricTensors.SymmetrizeTwoIndexTensor!(Mij2);
    SymmetricTensors.SymmetrizeThreeIndexTensor!(Mijk3);
    SymmetricTensors.SymmetrizeFourIndexTensor!(Mijkl4);
    SymmetricTensors.SymmetrizeTwoIndexTensor!(Sij2);
    SymmetricTensors.SymmetrizeThreeIndexTensor!(Sijk3);
    
    return t, Mij2, Mijk3, Mijkl4, Sij2, Sijk3
end

function load_trajectory(sol_filename::String)
    h5f = h5open(sol_filename, "r")
    t = h5f["t"][:]
    r = h5f["Trajectory/r"][:]
    θ = h5f["Trajectory/theta"][:]
    ϕ = h5f["Trajectory/phi"][:]
    # dr_dt = h5f["Trajectory/r_dot"][:]
    # dθ_dt = h5f["Trajectory/theta_dot"][:]
    # dϕ_dt = h5f["Trajectory/phi_dot"][:]
    # d2r_dt2 = h5f["Trajectory/r_ddot"][:]
    # d2θ_dt2 = h5f["Trajectory/theta_ddot"][:]
    # d2ϕ_dt2 = h5f["Trajectory/phi_ddot"][:]
    # dt_dτ = h5f["Trajectory/Gamma"][:]
    close(h5f)
    return t, r, θ, ϕ
end

function load_constants_of_motion(sol_filename::String)
    h5f = h5open(sol_filename, "r")
    t = h5f["ConstantsOfMotion/t"][:]
    EE = h5f["ConstantsOfMotion/Energy"][:]
    LL = h5f["ConstantsOfMotion/AngularMomentum"][:]
    CC = h5f["ConstantsOfMotion/CarterConstant"][:]
    QQ = h5f["ConstantsOfMotion/AltCarterConstant"][:]
    pArray = h5f["ConstantsOfMotion/p"][:]
    ecc = h5f["ConstantsOfMotion/eccentricity"][:]
    θminArray = h5f["ConstantsOfMotion/theta_min"][:]
    close(h5f)
    return t, EE, LL, QQ, CC, pArray, ecc, θminArray
end

function load_fluxes(sol_filename::String)
    h5f = h5open(sol_filename, "r")
    t = h5f["ConstantsOfMotion/t"][:]
    Edot = h5f["ConstantsOfMotion/Edot"][:]
    Ldot = h5f["ConstantsOfMotion/Ldot"][:]
    Qdot = h5f["ConstantsOfMotion/Qdot"][:]
    Cdot = h5f["ConstantsOfMotion/Cdot"][:]
    close(h5f)
    return t, Edot, Ldot, Qdot, Cdot
end

end