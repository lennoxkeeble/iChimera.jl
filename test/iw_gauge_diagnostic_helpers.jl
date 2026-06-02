using iChimera
using DifferentialEquations
using LinearAlgebra
using Printf
using StaticArrays
using Statistics

const SH = iChimera.SelfAccelerationHarmonic
const SHIW = iChimera.SelfAccelerationHarmonicIW

zero_mvector3() = @MVector [0.0 for _ in 1:3]
zero_mvector4() = @MVector [0.0 for _ in 1:4]
zero_mmatrix33() = @MMatrix [0.0 for _ in 1:3, _ in 1:3]
zero_marray333() = @MArray [0.0 for _ in 1:3, _ in 1:3, _ in 1:3]
zero_marray3333() = @MArray [0.0 for _ in 1:3, _ in 1:3, _ in 1:3, _ in 1:3]

function compute_bl_time_frequencies(a::Float64, p::Float64, e::Float64, θmin::Float64, sign_Lz::Int)
    ω = iChimera.ConstantsOfMotion.KerrFreqs(a, p, e, θmin, sign_Lz)
    return ω[1:3] ./ ω[4]
end

function compute_geodesic_duration(a::Float64, p::Float64, e::Float64, θmin::Float64, sign_Lz::Int, nOrbits::Float64)
    Ω = compute_bl_time_frequencies(a, p, e, θmin, sign_Lz)
    periods = Float64[]
    for Ωi in (Ω[1], Ω[2])
        if isfinite(Ωi) && Ωi != 0.0
            push!(periods, 2π / abs(Ωi))
        end
    end
    isempty(periods) && error("Could not define T_max from the BL-time radial/polar frequencies.")
    T_max = maximum(periods)
    return nOrbits * T_max, Ω, T_max
end

function allocate_flux_workspace()
    return (
        xBL = zero_mvector3(),
        vBL = zero_mvector3(),
        aBL = zero_mvector3(),
        dxBL_dt = zero_mvector3(),
        d2xBL_dt = zero_mvector3(),
        d3xBL_dt = zero_mvector3(),
        d4xBL_dt = zero_mvector3(),
        d5xBL_dt = zero_mvector3(),
        d6xBL_dt = zero_mvector3(),
        d7xBL_dt = zero_mvector3(),
        d8xBL_dt = zero_mvector3(),
        d9xBL_dt = zero_mvector3(),
        dx_dlambda = zero_mvector3(),
        d2x_dlambda = zero_mvector3(),
        d3x_dlambda = zero_mvector3(),
        d4x_dlambda = zero_mvector3(),
        d5x_dlambda = zero_mvector3(),
        d6x_dlambda = zero_mvector3(),
        d7x_dlambda = zero_mvector3(),
        d8x_dlambda = zero_mvector3(),
        d9x_dlambda = zero_mvector3(),
        xH = zero_mvector3(),
        dxH_dt = zero_mvector3(),
        d2xH_dt = zero_mvector3(),
        d3xH_dt = zero_mvector3(),
        d4xH_dt = zero_mvector3(),
        d5xH_dt = zero_mvector3(),
        d6xH_dt = zero_mvector3(),
        d7xH_dt = zero_mvector3(),
        d8xH_dt = zero_mvector3(),
        d9xH_dt = zero_mvector3(),
        vH = zero_mvector3(),
        aH = zero_mvector3(),
        aSF_H_old = zero_mvector4(),
        aSF_BL_old = zero_mvector4(),
        aSF_H_new = zero_mvector4(),
        aSF_BL_new = zero_mvector4(),
        Virr = zero_mvector3(),
        dVrr_da = zero_mvector3(),
        dVirr_dt = zero_mvector3(),
        dVirr_da = zero_mmatrix33(),
        Mij5 = zero_mmatrix33(),
        Mij6 = zero_mmatrix33(),
        Mij7 = zero_mmatrix33(),
        Mij8 = zero_mmatrix33(),
        dxmMij5 = zero_marray333(),
        dxmMij6 = zero_marray333(),
        dxmMij7 = zero_marray333(),
        Mijk7 = zero_marray333(),
        Mijk8 = zero_marray333(),
        dxmMijk7 = zero_marray3333(),
        Sij5 = zero_mmatrix33(),
        Sij6 = zero_mmatrix33(),
        dxmSij5 = zero_marray333(),
    )
end

function zero_higher_multipoles!(ws)
    ws.Mijk7 .= 0.0
    ws.Mijk8 .= 0.0
    ws.dxmMijk7 .= 0.0
    ws.Sij5 .= 0.0
    ws.Sij6 .= 0.0
    ws.dxmSij5 .= 0.0
    return nothing
end

function build_geodesic_integrator(;
    a::Float64,
    p::Float64,
    e::Float64,
    θmin::Float64,
    sign_Lz::Int,
    psi0::Float64=0.0,
    chi0::Float64=0.0,
    phi0::Float64=0.0,
    nOrbits::Float64=1.0,
    reltol::Float64=1.0e-13,
    abstol::Float64=1.0e-13,
    maxiters::Int=Int(1e8),
)
    t_geodesic, Ω, T_max = compute_geodesic_duration(a, p, e, θmin, sign_Lz, nOrbits)
    E, L, Q, C, _, p3, p4, zp, zm =
        iChimera.BLTimeGeodesics.compute_ODE_params(a, p, e, θmin, sign_Lz)
    params = @SArray [a, E, L, p, e, θmin, p3, p4, zp, zm]
    ics = @SArray [psi0, chi0, phi0]
    tspan = (0.0, t_geodesic)
    prob = e == 0.0 ?
        ODEProblem(iChimera.BLTimeGeodesics.HJ_Eqns_circular, ics, tspan, params) :
        ODEProblem(iChimera.BLTimeGeodesics.HJ_Eqns, ics, tspan, params)
    dt0 = t_geodesic / 1000.0
    integrator = init(
        prob,
        AutoTsit5(RK4()),
        adaptive=true,
        dt=dt0,
        reltol=reltol,
        abstol=abstol,
        maxiters=maxiters,
    )
    return (
        integrator=integrator,
        E=E,
        L=L,
        Q=Q,
        C=C,
        p3=p3,
        p4=p4,
        zp=zp,
        zm=zm,
        Ω=Ω,
        T_max=T_max,
        t_geodesic=t_geodesic,
    )
end

function step_integrator_to!(integrator, target_t::Float64)
    if target_t > integrator.t
        step!(integrator, target_t - integrator.t, true)
    end
    return integrator
end

function build_local_state!(
    ws,
    integrator,
    a::Float64,
    p::Float64,
    e::Float64,
    θmin::Float64,
    q::Float64,
    E::Float64,
    L::Float64,
    Q::Float64,
    C::Float64,
    p3::Float64,
    p4::Float64,
    zp::Float64,
    zm::Float64;
    onepn::Float64=0.0,
    twopn::Float64=0.0,
    twop5pn::Float64=0.0,
    v2power::Int=0,
)
    tt, rr, θθ, ϕϕ, r_dot, θ_dot, ϕ_dot, r_ddot, θ_ddot, ϕ_ddot, dt_dτ, _, _ =
        iChimera.Inspiral.compute_geodesic_arrays(integrator, a, E, L, p, e, θmin, p3, p4, zp, zm)

    ws.xBL[1] = rr
    ws.xBL[2] = θθ
    ws.xBL[3] = ϕϕ
    ws.vBL[1] = r_dot
    ws.vBL[2] = θ_dot
    ws.vBL[3] = ϕ_dot
    ws.aBL[1] = r_ddot
    ws.aBL[2] = θ_ddot
    ws.aBL[3] = ϕ_ddot

    iChimera.CoordinateDerivs.ComputeDerivs!(
        ws.xBL,
        sign(ws.vBL[1]),
        sign(ws.vBL[2]),
        ws.dxBL_dt,
        ws.d2xBL_dt,
        ws.d3xBL_dt,
        ws.d4xBL_dt,
        ws.d5xBL_dt,
        ws.d6xBL_dt,
        ws.d7xBL_dt,
        ws.d8xBL_dt,
        ws.d9xBL_dt,
        ws.dx_dlambda,
        ws.d2x_dlambda,
        ws.d3x_dlambda,
        ws.d4x_dlambda,
        ws.d5x_dlambda,
        ws.d6x_dlambda,
        ws.d7x_dlambda,
        ws.d8x_dlambda,
        ws.d9x_dlambda,
        a,
        E,
        L,
        C,
    )

    iChimera.HarmonicCoordDerivs.compute_harmonic_derivs!(
        ws.xBL,
        ws.dxBL_dt,
        ws.d2xBL_dt,
        ws.d3xBL_dt,
        ws.d4xBL_dt,
        ws.d5xBL_dt,
        ws.d6xBL_dt,
        ws.d7xBL_dt,
        ws.d8xBL_dt,
        ws.d9xBL_dt,
        ws.xH,
        ws.dxH_dt,
        ws.d2xH_dt,
        ws.d3xH_dt,
        ws.d4xH_dt,
        ws.d5xH_dt,
        ws.d6xH_dt,
        ws.d7xH_dt,
        ws.d8xH_dt,
        ws.d9xH_dt,
        a,
    )

    iChimera.MultipoleDerivs.compute_SF_moments!(
        q,
        ws.Mij5,
        ws.Mij6,
        ws.Mij7,
        ws.Mij8,
        ws.dxmMij5,
        ws.dxmMij6,
        ws.dxmMij7,
        ws.Mijk7,
        ws.Mijk8,
        ws.dxmMijk7,
        ws.Sij5,
        ws.Sij6,
        ws.dxmSij5,
        ws.xH,
        ws.dxH_dt,
        ws.d2xH_dt,
        ws.d3xH_dt,
        ws.d4xH_dt,
        ws.d5xH_dt,
        ws.d6xH_dt,
        ws.d7xH_dt,
        ws.d8xH_dt,
        ws.d9xH_dt,
        onepn,
        twopn,
        twop5pn,
    )

    Vrr, dVrr_dt = iChimera.RRPotentials.compute_RR_potentials!(
        ws.Virr,
        ws.dVrr_da,
        ws.dVirr_dt,
        ws.dVirr_da,
        ws.xH,
        ws.dxH_dt,
        ws.d2xH_dt,
        ws.Mij5,
        ws.Mij6,
        ws.Mij7,
        ws.Mij8,
        ws.dxmMij5,
        ws.dxmMij6,
        ws.dxmMij7,
        ws.Mijk7,
        ws.Mijk8,
        ws.dxmMijk7,
        ws.Sij5,
        ws.Sij6,
        ws.dxmSij5,
        v2power,
    )

    ws.vH .= ws.dxH_dt
    ws.aH .= ws.d2xH_dt
    rH = SH.norm_3d(ws.xH)
    v = SH.norm_3d(ws.vH)

    return (
        tt=tt,
        rr=rr,
        θθ=θθ,
        ϕϕ=ϕϕ,
        r_dot=r_dot,
        θ_dot=θ_dot,
        ϕ_dot=ϕ_dot,
        dt_dτ=dt_dτ,
        Vrr=Vrr,
        dVrr_dt=dVrr_dt,
        rH=rH,
        v=v,
    )
end

function compute_self_force_pair!(
    ws,
    local_state,
    spin::Float64,
    q::Float64;
    alpha_current::Float64=4.0,
    beta_current::Float64=5.0,
    alpha_target::Float64=4.0,
    beta_target::Float64=5.0,
    gauge_sign::Float64=-1.0,
    eps_fd::Float64=1.0e-5 * max(1.0, SH.norm_3d(ws.xH)),
)
    SH.aRRα(
        ws.aSF_H_old,
        ws.aSF_BL_old,
        ws.xH,
        local_state.v,
        ws.vH,
        ws.xBL,
        local_state.rH,
        spin,
        local_state.Vrr,
        local_state.dVrr_dt,
        ws.Virr,
        ws.dVrr_da,
        ws.dVirr_dt,
        ws.dVirr_da,
    )

    SHIW.aRRα(
        ws.aSF_H_new,
        ws.aSF_BL_new,
        ws.xH,
        local_state.v,
        ws.vH,
        ws.aH,
        ws.d3xH_dt,
        ws.xBL,
        local_state.rH,
        spin,
        q,
        local_state.Vrr,
        local_state.dVrr_dt,
        ws.Virr,
        ws.dVrr_da,
        ws.dVirr_dt,
        ws.dVirr_da;
        alpha_current=alpha_current,
        beta_current=beta_current,
        alpha_target=alpha_target,
        beta_target=beta_target,
        gauge_sign=gauge_sign,
        eps_fd=eps_fd,
    )

    return nothing
end

function flux_triplet(
    spin::Float64,
    rr::Float64,
    θθ::Float64,
    ϕϕ::Float64,
    dt_dτ::Float64,
    r_dot::Float64,
    θ_dot::Float64,
    ϕ_dot::Float64,
    aSF_BL,
    E::Float64,
    L::Float64,
    Q::Float64,
    C::Float64,
    p::Float64,
    e::Float64,
    θmin::Float64,
)
    _, Edot, _, Ldot, _, Qdot, _, Cdot, _, _, _ = iChimera.EvolveConstants.Evolve_BL(
        0.0,
        spin,
        rr,
        θθ,
        ϕϕ,
        dt_dτ,
        r_dot,
        θ_dot,
        ϕ_dot,
        aSF_BL,
        E,
        0.0,
        L,
        0.0,
        Q,
        0.0,
        C,
        0.0,
        p,
        e,
        θmin,
    )
    return (Edot=Edot, Ldot=Ldot, Qdot=Qdot, Cdot=Cdot)
end

function compute_local_flux_pair!(
    ws,
    integrator,
    spin::Float64,
    p::Float64,
    e::Float64,
    θmin::Float64,
    q::Float64,
    E::Float64,
    L::Float64,
    Q::Float64,
    C::Float64,
    p3::Float64,
    p4::Float64,
    zp::Float64,
    zm::Float64;
    onepn::Float64=0.0,
    twopn::Float64=0.0,
    twop5pn::Float64=0.0,
    v2power::Int=0,
    alpha_current::Float64=4.0,
    beta_current::Float64=5.0,
    alpha_target::Float64=4.0,
    beta_target::Float64=5.0,
    gauge_sign::Float64=-1.0,
    eps_fd::Union{Nothing, Float64}=nothing,
)
    local_state = build_local_state!(
        ws,
        integrator,
        spin,
        p,
        e,
        θmin,
        q,
        E,
        L,
        Q,
        C,
        p3,
        p4,
        zp,
        zm;
        onepn=onepn,
        twopn=twopn,
        twop5pn=twop5pn,
        v2power=v2power,
    )

    eps_use = isnothing(eps_fd) ? 1.0e-5 * max(1.0, SH.norm_3d(ws.xH)) : eps_fd
    compute_self_force_pair!(
        ws,
        local_state,
        spin,
        q;
        alpha_current=alpha_current,
        beta_current=beta_current,
        alpha_target=alpha_target,
        beta_target=beta_target,
        gauge_sign=gauge_sign,
        eps_fd=eps_use,
    )

    old_flux = flux_triplet(
        spin,
        local_state.rr,
        local_state.θθ,
        local_state.ϕϕ,
        local_state.dt_dτ,
        local_state.r_dot,
        local_state.θ_dot,
        local_state.ϕ_dot,
        ws.aSF_BL_old,
        E,
        L,
        Q,
        C,
        p,
        e,
        θmin,
    )
    new_flux = flux_triplet(
        spin,
        local_state.rr,
        local_state.θθ,
        local_state.ϕϕ,
        local_state.dt_dτ,
        local_state.r_dot,
        local_state.θ_dot,
        local_state.ϕ_dot,
        ws.aSF_BL_new,
        E,
        L,
        Q,
        C,
        p,
        e,
        θmin,
    )

    return local_state, old_flux, new_flux
end

function abs_rel_diff(x, y)
    diff = norm(x - y)
    scale = max(norm(x), eps(Float64))
    return diff, diff / scale
end

function dummy_rr_data(; scale::Float64=1.0)
    Vrr = 1.7e-9 * scale
    dVrr_dt = -2.3e-11 * scale
    Virr = @MVector [4.0e-10 * scale, -7.0e-10 * scale, 9.0e-10 * scale]
    dVrr_da = @MVector [3.0e-11 * scale, -4.0e-11 * scale, 8.0e-11 * scale]
    dVirr_dt = @MVector [5.0e-12 * scale, -6.0e-12 * scale, 7.0e-12 * scale]
    dVirr_da = zero_mmatrix33()
    dVirr_da[1, 1] = 9.0e-13 * scale
    dVirr_da[1, 2] = -2.0e-12 * scale
    dVirr_da[1, 3] = 3.0e-12 * scale
    dVirr_da[2, 1] = 4.0e-12 * scale
    dVirr_da[2, 2] = 5.0e-13 * scale
    dVirr_da[2, 3] = -6.0e-12 * scale
    dVirr_da[3, 1] = -7.0e-12 * scale
    dVirr_da[3, 2] = 8.0e-12 * scale
    dVirr_da[3, 3] = 9.0e-13 * scale

    return (; Vrr, dVrr_dt, Virr, dVrr_da, dVirr_dt, dVirr_da)
end

function zero_rr_data()
    Vrr = 0.0
    dVrr_dt = 0.0
    Virr = zero_mvector3()
    dVrr_da = zero_mvector3()
    dVirr_dt = zero_mvector3()
    dVirr_da = zero_mmatrix33()
    return (; Vrr, dVrr_dt, Virr, dVrr_da, dVirr_dt, dVirr_da)
end

function newtonian_jerk(xH, vH)
    R = norm(xH)
    return -vH / R^3 + 3.0 * xH * dot(xH, vH) / R^5
end

function expected_IW_shift(xH, vH, q, alpha_current, beta_current, alpha_target, beta_target)
    R = norm(xH)
    n = xH ./ R
    vr = dot(n, vH)
    v2 = dot(vH, vH)
    Δα = alpha_target - alpha_current
    Δβ = beta_target - beta_current
    ΔA = Δβ * (3.0 * v2 - 3.0 / R - 5.0 * vr^2) + 2.0 * Δα / R
    ΔB = Δα * (v2 - 1.0 / R - 3.0 * vr^2)
    return (8.0 / 5.0) * q / R^3 .* (ΔA * vr .* n .- ΔB .* vH)
end

function sf_to_coordinate_accel(aSF_H, vH, xH, spin)
    ΓH = SH.Γ(vH, xH, spin)
    return [(aSF_H[i + 1] - vH[i] * aSF_H[1]) / ΓH^2 for i in 1:3]
end

function self_force_difference_to_coord_accel(aSF_H_new, aSF_H_old, vH, xH, spin)
    δaSF = aSF_H_new .- aSF_H_old
    return sf_to_coordinate_accel(δaSF, vH, xH, spin)
end

function print_repo_summary()
    println("Repository gauge/flux call-path summary:")
    println("  potentials: src/RRPotentials.jl :: compute_RR_potentials!")
    println("  harmonic derivatives: src/HarmonicCoordDerivs.jl :: compute_harmonic_derivs!")
    println("  old harmonic self-force: src/SelfAccelerationHarmonic.jl :: aRRα")
    println("  IW harmonic self-force: src/SelfAccelerationHarmonicIW.jl :: aRRα")
    println("  harmonic metric map: src/HarmonicCoords.jl :: g_μν_H / gμν_H")
    println("  BL->harmonic acceleration map: src/HarmonicCoords.jl :: aHtoBL")
    println("  BL-time frequencies: src/ConstantsOfMotion.jl :: KerrFreqs")
    println("  fixed-geodesic harmonic call site: src/Inspiral.jl and test/compare_old_new_fluxes.jl")
    return nothing
end
