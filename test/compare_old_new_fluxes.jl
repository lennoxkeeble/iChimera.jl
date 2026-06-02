using iChimera
using DifferentialEquations
using StaticArrays

module PlotSettings
using CairoMakie

function set_global_themes!()
    set_theme!(theme_latexfonts())
    return nothing
end

load_settings() = (900, 600, 1400, 600, 900, 900, 28, 30, 30, 24, 24, 25, true, true)

end

include(joinpath(@__DIR__, "..", "src", "QuickPlots.jl"))

zero_mvector3() = @MVector [0.0 for i in 1:3]
zero_mvector4() = @MVector [0.0 for i in 1:4]
zero_mmatrix33() = @MMatrix [0.0 for i in 1:3, j in 1:3]
zero_marray333() = @MArray [0.0 for i in 1:3, j in 1:3, k in 1:3]
zero_marray3333() = @MArray [0.0 for i in 1:3, j in 1:3, k in 1:3, l in 1:3]

function compute_bl_time_frequencies(a::Float64, p::Float64, e::Float64, θmin::Float64, sign_Lz::Int)
    ω = iChimera.ConstantsOfMotion.KerrFreqs(a, p, e, θmin, sign_Lz)
    return ω[1:3] ./ ω[4]
end

function compute_geodesic_duration(a::Float64, p::Float64, e::Float64, θmin::Float64, sign_Lz::Int, nOrbits::Float64)
    Ω = compute_bl_time_frequencies(a, p, e, θmin, sign_Lz)

    periods = Float64[]
    for Ωi in (Ω[1], Ω[2])
        if abs(Ωi) < 1e9 && Ωi != 0.0
            push!(periods, 2π / abs(Ωi))
        end
    end

    isempty(periods) && error("Could not define T_max from finite BL-time radial/polar frequencies.")

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
        dxmSij5 = zero_marray333()
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

function compute_local_edot_pair!(
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
    alpha_current::Float64=4.0,
    beta_current::Float64=5.0,
    alpha_target::Float64=4.0,
    beta_target::Float64=5.0,
    gauge_sign::Float64=-1.0,
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

    # Mass quadrupole only, with all PN corrections disabled.
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
        0.0,
        0.0,
        0.0,
    )
    zero_higher_multipoles!(ws)

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
        0,
    )

    ws.vH .= ws.dxH_dt
    ws.aH .= ws.d2xH_dt
    rH = iChimera.SelfAccelerationHarmonic.norm_3d(ws.xH)
    v = iChimera.SelfAccelerationHarmonic.norm_3d(ws.vH)

    iChimera.SelfAccelerationHarmonic.aRRα(
        ws.aSF_H_old,
        ws.aSF_BL_old,
        ws.xH,
        v,
        ws.vH,
        ws.xBL,
        rH,
        a,
        Vrr,
        dVrr_dt,
        ws.Virr,
        ws.dVrr_da,
        ws.dVirr_dt,
        ws.dVirr_da,
    )

    iChimera.SelfAccelerationHarmonicIW.aRRα(
        ws.aSF_H_new,
        ws.aSF_BL_new,
        ws.xH,
        v,
        ws.vH,
        ws.aH,
        ws.d3xH_dt,
        ws.xBL,
        rH,
        a,
        q,
        Vrr,
        dVrr_dt,
        ws.Virr,
        ws.dVrr_da,
        ws.dVirr_dt,
        ws.dVirr_da;
        alpha_current=alpha_current,
        beta_current=beta_current,
        alpha_target=alpha_target,
        beta_target=beta_target,
        gauge_sign=gauge_sign,
        eps_fd=1.0e-5 * max(1.0, iChimera.SelfAccelerationHarmonic.norm_3d(ws.xH)),
    )

    _, Edot_old, _, _, _, _, _, _, _, _, _ = iChimera.EvolveConstants.Evolve_BL(
        0.0,
        a,
        rr,
        θθ,
        ϕϕ,
        dt_dτ,
        r_dot,
        θ_dot,
        ϕ_dot,
        ws.aSF_BL_old,
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

    _, Edot_new, _, _, _, _, _, _, _, _, _ = iChimera.EvolveConstants.Evolve_BL(
        0.0,
        a,
        rr,
        θθ,
        ϕϕ,
        dt_dτ,
        r_dot,
        θ_dot,
        ϕ_dot,
        ws.aSF_BL_new,
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

    return tt, Edot_old, Edot_new
end

function compare_old_new_fluxes(;
    a::Float64=0.8,
    p::Float64=10.0,
    e::Float64=0.5,
    θmin::Float64=π / 3,
    sign_Lz::Int=1,
    q::Float64=1.0e-5,
    psi0::Float64=0.0,
    chi0::Float64=0.0,
    phi0::Float64=0.0,
    nOrbits::Float64=4.0,
    n_samples::Int=400,
    reltol::Float64=1.0e-13,
    abstol::Float64=1.0e-13,
    maxiters::Int=Int(1e8),
    alpha_current::Float64=4.0,
    beta_current::Float64=5.0,
    alpha_target::Float64=4.0,
    beta_target::Float64=5.0,
    gauge_sign::Float64=-1.0,
)
    n_samples < 2 && error("n_samples must be at least 2.")
    nOrbits <= 0.0 && error("nOrbits must be positive.")
    iw_active = !(alpha_current == alpha_target && beta_current == beta_target)

    t_geodesic, Ω, T_max = compute_geodesic_duration(a, p, e, θmin, sign_Lz, nOrbits)

    E, L, Q, C, _, p3, p4, zp, zm =
        iChimera.BLTimeGeodesics.compute_ODE_params(a, p, e, θmin, sign_Lz)

    params = @SArray [a, E, L, p, e, θmin, p3, p4, zp, zm]
    ics = @SArray [psi0, chi0, phi0]
    tspan = (0.0, t_geodesic)

    prob = e == 0.0 ?
        ODEProblem(iChimera.BLTimeGeodesics.HJ_Eqns_circular, ics, tspan, params) :
        ODEProblem(iChimera.BLTimeGeodesics.HJ_Eqns, ics, tspan, params)

    sample_times = collect(range(0.0, t_geodesic, length=n_samples))
    initial_dt = (sample_times[2] - sample_times[1]) / 1000.0
    integrator = init(
        prob,
        AutoTsit5(RK4()),
        adaptive=true,
        dt=initial_dt,
        reltol=reltol,
        abstol=abstol,
        maxiters=maxiters,
    )

    ws = allocate_flux_workspace()
    Edot_old = zeros(n_samples)
    Edot_new = zeros(n_samples)

    for (idx, t_sample) in enumerate(sample_times)
        if idx > 1
            step!(integrator, t_sample - integrator.t, true)
        end

        _, Edot_old[idx], Edot_new[idx] = compute_local_edot_pair!(
            ws,
            integrator,
            a,
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
            alpha_current=alpha_current,
            beta_current=beta_current,
            alpha_target=alpha_target,
            beta_target=beta_target,
            gauge_sign=gauge_sign,
        )
    end

    return (
        times=sample_times,
        Edot_old=Edot_old,
        Edot_new=Edot_new,
        Ω=Ω,
        T_max=T_max,
        t_geodesic=t_geodesic,
        nOrbits=nOrbits,
        iw_active=iw_active,
        alpha_current=alpha_current,
        beta_current=beta_current,
        alpha_target=alpha_target,
        beta_target=beta_target,
        gauge_sign=gauge_sign,
    )
end

# By default the IW wrapper is evaluated in the same Burke-Thorne gauge as the
# legacy harmonic path, so the two curves should lie on top of one another. To
# test a nontrivial gauge mapping, change alpha_current/beta_current while
# keeping alpha_target=4 and beta_target=5.

# Interactive comparison settings.
# The legacy harmonic curve always uses the existing SelfAccelerationHarmonic
# infrastructure. The IW curve below uses the gauged SelfAccelerationHarmonicIW
# path with the explicit current -> target gauge map chosen here.
alpha_current_compare = 3.0
beta_current_compare = 4.0
alpha_target_compare = 4.0
beta_target_compare = 5.0
gauge_sign_compare = -1.0

flux_comparison = compare_old_new_fluxes(
    alpha_current=alpha_current_compare,
    beta_current=beta_current_compare,
    alpha_target=alpha_target_compare,
    beta_target=beta_target_compare,
    gauge_sign=gauge_sign_compare,
)

lim_x_min = nothing;
lim_x_max = nothing;
lim_y_min = nothing;
lim_y_max = nothing;

QuickPlots.plot11(
    [flux_comparison.times, flux_comparison.times],
    [flux_comparison.Edot_old, flux_comparison.Edot_new];
    colors=[:royalblue, :tomato],
    labels=[
        "old harmonic",
        "IW gauge map ($(flux_comparison.alpha_current), $(flux_comparison.beta_current)) → ($(flux_comparison.alpha_target), $(flux_comparison.beta_target))",
    ],
    linestyles=[:solid, :dash],
    linewidths=[3.0, 3.0],
    xlabel=L"t / M",
    ylabel=L"\dot{E}",
    text="mass quadrupole only\nIW active = $(flux_comparison.iw_active)\n(α,β): ($(flux_comparison.alpha_current), $(flux_comparison.beta_current)) → ($(flux_comparison.alpha_target), $(flux_comparison.beta_target))\nnOrbits = $(flux_comparison.nOrbits)\nT_max = $(round(flux_comparison.T_max; digits=6)) M",
    text_xloc=0.05,
    text_yloc=0.95,
    legend=true,
    position=:rb,
    lim_x_min = lim_x_min,
    lim_x_max = lim_x_max,
    lim_y_min = lim_y_min,
    lim_y_max = lim_y_max,
)
