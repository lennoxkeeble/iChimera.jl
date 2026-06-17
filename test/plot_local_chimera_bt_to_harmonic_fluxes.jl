using iChimera
using DifferentialEquations
using StaticArrays
using LinearAlgebra
using LaTeXStrings

module PlotSettings
using CairoMakie

function set_global_themes!()
    set_theme!(theme_latexfonts())
    return nothing
end

load_settings() = (900, 600, 1400, 600, 900, 900, 28, 30, 30, 24, 24, 25, true, true)

end

include(joinpath(@__DIR__, "..", "src", "QuickPlots.jl"))

const HRR = iChimera.HarmonicRRAcceleration

"""
User-editable parameters.

This script samples the local fluxes on a fixed Kerr geodesic using the hybrid
BT->harmonic correction model:

    Chimera + Δa_BT->H

It compares two versions of the correction:

1. Chimera + 2.5PN BT->harmonic correction only
2. Chimera + 2.5PN + 3.5PN BT->harmonic correction

The geodesic duration is

    t_span = time_fraction_of_max_period * max(T_r, T_θ),

with periods defined using BL coordinate-time frequencies.
"""

# ----------------------------- EMRI / geodesic parameters -----------------------------
a = 0.7
p = 10.0
e = 0.2
q = 1.0e-5

inclination_deg = 45.0
inclination_type = "iota"          # "iota" or "theta_inc"
sign_Lz = inclination_deg < 90.0 ? 1 : -1

use_direct_theta_min = false
theta_min_direct = π / 3           # used only if use_direct_theta_min = true

psi0 = 0.1
chi0 = 0.2
phi0 = 0.3

lmax_mass_fluxes = 3
lmax_current_fluxes = 2
OnePN = 1.0
TwoPN = 1.0
TwoPointFivePN = 1.0

time_fraction_of_max_period = 10.0
n_samples = 400

reltol = 1.0e-13
abstol = 1.0e-13
maxiters = Int(1e8)

save_plots = true
results_dir = joinpath(@__DIR__, "..", "Results", "local_chimera_bt_to_harmonic_fluxes")
plot_tag = "a$(a)_p$(p)_e$(e)_q$(q)"
# --------------------------------------------------------------------------------------

function compute_theta_min_for_script(a::Float64, p::Float64, e::Float64, inclination_deg::Float64, inclination_type::String, sign_Lz::Int64;
                                      use_direct_theta_min::Bool=false, theta_min_direct::Float64=π / 2)
    if use_direct_theta_min
        return theta_min_direct
    end
    return iChimera.compute_theta_min(a, p, e, inclination_deg, inclination_type, sign_Lz)
end

function compute_bl_time_frequencies(a::Float64, p::Float64, e::Float64, θmin::Float64, sign_Lz::Int64)
    ω = iChimera.ConstantsOfMotion.KerrFreqs(a, p, e, θmin, sign_Lz)
    return ω[1:3] ./ ω[4]
end

function compute_geodesic_duration(a::Float64, p::Float64, e::Float64, θmin::Float64, sign_Lz::Int64, time_fraction_of_max_period::Float64)
    time_fraction_of_max_period > 0.0 || error("time_fraction_of_max_period must be positive.")

    Ω = compute_bl_time_frequencies(a, p, e, θmin, sign_Lz)
    periods = Float64[]

    for Ωi in (Ω[1], Ω[2])
        if Ωi != 0.0 && abs(Ωi) < 1.0e9
            push!(periods, 2π / abs(Ωi))
        end
    end

    isempty(periods) && error("Could not define a finite BL-time max period from Ω_r and Ω_θ. Avoid a fully degenerate circular-equatorial configuration for this script.")

    T_max = maximum(periods)
    return time_fraction_of_max_period * T_max, Ω, T_max
end

function local_fluxes_from_bl_four_accel(a::Float64, p::Float64, e::Float64, θmin::Float64,
                                         E::Float64, L::Float64, Q::Float64, C::Float64,
                                         rr::Float64, θθ::Float64, ϕϕ::Float64,
                                         dt_dτ::Float64, r_dot::Float64, θ_dot::Float64, ϕ_dot::Float64,
                                         aSF_BL::AbstractVector)
    _, dE_dt, _, dL_dt, _, dQ_dt, _, _, _, _, _ = iChimera.EvolveConstants.Evolve_BL(
        0.0,
        a,
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

    return dE_dt, dL_dt, dQ_dt
end

function allocate_rr_workspace()
    return (
        xBL=@MArray zeros(3),
        vBL=@MArray zeros(3),
        aBL=@MArray zeros(3),
        dxBL_dt=@MArray zeros(3),
        d2xBL_dt=@MArray zeros(3),
        d3xBL_dt=@MArray zeros(3),
        d4xBL_dt=@MArray zeros(3),
        d5xBL_dt=@MArray zeros(3),
        d6xBL_dt=@MArray zeros(3),
        d7xBL_dt=@MArray zeros(3),
        d8xBL_dt=@MArray zeros(3),
        d9xBL_dt=@MArray zeros(3),
        dx_dλ=@MArray zeros(3),
        d2x_dλ=@MArray zeros(3),
        d3x_dλ=@MArray zeros(3),
        d4x_dλ=@MArray zeros(3),
        d5x_dλ=@MArray zeros(3),
        d6x_dλ=@MArray zeros(3),
        d7x_dλ=@MArray zeros(3),
        d8x_dλ=@MArray zeros(3),
        d9x_dλ=@MArray zeros(3),
        xH=@MArray zeros(3),
        dxH_dt=@MArray zeros(3),
        d2xH_dt=@MArray zeros(3),
        d3xH_dt=@MArray zeros(3),
        d4xH_dt=@MArray zeros(3),
        d5xH_dt=@MArray zeros(3),
        d6xH_dt=@MArray zeros(3),
        d7xH_dt=@MArray zeros(3),
        d8xH_dt=@MArray zeros(3),
        d9xH_dt=@MArray zeros(3),
        Mij5=@MArray zeros(3, 3),
        Mij6=@MArray zeros(3, 3),
        Mij7=@MArray zeros(3, 3),
        Mij8=@MArray zeros(3, 3),
        dxmMij5=@MArray zeros(3, 3, 3),
        dxmMij6=@MArray zeros(3, 3, 3),
        dxmMij7=@MArray zeros(3, 3, 3),
        Mijk7=@MArray zeros(3, 3, 3),
        Mijk8=@MArray zeros(3, 3, 3),
        dxmMijk7=@MArray zeros(3, 3, 3, 3),
        Sij5=@MArray zeros(3, 3),
        Sij6=@MArray zeros(3, 3),
        dxmSij5=@MArray zeros(3, 3, 3),
        aSF_H=@MArray zeros(4),
        aSF_BL=@MArray zeros(4),
        aSF25_H=@MArray zeros(4),
        aSF25_BL=@MArray zeros(4),
        aSF35_H=@MArray zeros(4),
        aSF35_BL=@MArray zeros(4),
        aSF_corr_H=@MArray zeros(4),
        aSF_corr_BL=@MArray zeros(4),
        Acoord_BT=@MArray zeros(3),
        Acoord_corr_H=@MArray zeros(3),
        Virr=@MArray zeros(3),
        ∂Vrr_∂a=@MArray zeros(3),
        ∂Virr_∂t=@MArray zeros(3),
        ∂Virr_∂a=@MArray zeros(3, 3),
    )
end

function copy_four_vector!(dest::AbstractVector, src::AbstractVector)
    @inbounds for i in eachindex(dest, src)
        dest[i] = src[i]
    end
    return dest
end

function compute_local_flux_pair!(workspace,
                                  integrator,
                                  a::Float64, p::Float64, e::Float64, θmin::Float64, q::Float64,
                                  E::Float64, L::Float64, Q::Float64, C::Float64,
                                  p3::Float64, p4::Float64, zp::Float64, zm::Float64;
                                  lmax_mass_fluxes::Int64,
                                  lmax_current_fluxes::Int64,
                                  OnePN::Float64,
                                  TwoPN::Float64,
                                  TwoPointFivePN::Float64)
    tt, rr, θθ, ϕϕ, r_dot, θ_dot, ϕ_dot, r_ddot, θ_ddot, ϕ_ddot, dt_dτ, _, _ =
        iChimera.Inspiral.compute_geodesic_arrays(integrator, a, E, L, p, e, θmin, p3, p4, zp, zm)

    workspace.xBL[1] = rr
    workspace.xBL[2] = θθ
    workspace.xBL[3] = ϕϕ
    workspace.vBL[1] = r_dot
    workspace.vBL[2] = θ_dot
    workspace.vBL[3] = ϕ_dot
    workspace.aBL[1] = r_ddot
    workspace.aBL[2] = θ_ddot
    workspace.aBL[3] = ϕ_ddot

    iChimera.CoordinateDerivs.ComputeDerivs!(
        workspace.xBL,
        sign(workspace.vBL[1]),
        sign(workspace.vBL[2]),
        workspace.dxBL_dt,
        workspace.d2xBL_dt,
        workspace.d3xBL_dt,
        workspace.d4xBL_dt,
        workspace.d5xBL_dt,
        workspace.d6xBL_dt,
        workspace.d7xBL_dt,
        workspace.d8xBL_dt,
        workspace.d9xBL_dt,
        workspace.dx_dλ,
        workspace.d2x_dλ,
        workspace.d3x_dλ,
        workspace.d4x_dλ,
        workspace.d5x_dλ,
        workspace.d6x_dλ,
        workspace.d7x_dλ,
        workspace.d8x_dλ,
        workspace.d9x_dλ,
        a,
        E,
        L,
        C,
    )

    iChimera.HarmonicCoordDerivs.compute_harmonic_derivs!(
        workspace.xBL,
        workspace.dxBL_dt,
        workspace.d2xBL_dt,
        workspace.d3xBL_dt,
        workspace.d4xBL_dt,
        workspace.d5xBL_dt,
        workspace.d6xBL_dt,
        workspace.d7xBL_dt,
        workspace.d8xBL_dt,
        workspace.d9xBL_dt,
        workspace.xH,
        workspace.dxH_dt,
        workspace.d2xH_dt,
        workspace.d3xH_dt,
        workspace.d4xH_dt,
        workspace.d5xH_dt,
        workspace.d6xH_dt,
        workspace.d7xH_dt,
        workspace.d8xH_dt,
        workspace.d9xH_dt,
        a,
    )

    iChimera.MultipoleDerivs.compute_SF_moments!(
        q,
        workspace.Mij5,
        workspace.Mij6,
        workspace.Mij7,
        workspace.Mij8,
        workspace.dxmMij5,
        workspace.dxmMij6,
        workspace.dxmMij7,
        workspace.Mijk7,
        workspace.Mijk8,
        workspace.dxmMijk7,
        workspace.Sij5,
        workspace.Sij6,
        workspace.dxmSij5,
        workspace.xH,
        workspace.dxH_dt,
        workspace.d2xH_dt,
        workspace.d3xH_dt,
        workspace.d4xH_dt,
        workspace.d5xH_dt,
        workspace.d6xH_dt,
        workspace.d7xH_dt,
        workspace.d8xH_dt,
        workspace.d9xH_dt,
        OnePN,
        TwoPN,
        TwoPointFivePN,
    )

    if lmax_mass_fluxes != 3
        workspace.Mijk7 .= 0.0
        workspace.Mijk8 .= 0.0
        workspace.dxmMijk7 .= 0.0
    end

    if lmax_current_fluxes != 2
        workspace.Sij5 .= 0.0
        workspace.Sij6 .= 0.0
        workspace.dxmSij5 .= 0.0
    end

    v2power = 0
    Vrr, ∂Vrr_∂t = iChimera.RRPotentials.compute_RR_potentials!(
        workspace.Virr,
        workspace.∂Vrr_∂a,
        workspace.∂Virr_∂t,
        workspace.∂Virr_∂a,
        workspace.xH,
        workspace.dxH_dt,
        workspace.d2xH_dt,
        workspace.Mij5,
        workspace.Mij6,
        workspace.Mij7,
        workspace.Mij8,
        workspace.dxmMij5,
        workspace.dxmMij6,
        workspace.dxmMij7,
        workspace.Mijk7,
        workspace.Mijk8,
        workspace.dxmMijk7,
        workspace.Sij5,
        workspace.Sij6,
        workspace.dxmSij5,
        v2power,
    )

    iChimera.Inspiral.compute_chimera_self_accel!(
        workspace.aSF_H,
        workspace.aSF_BL,
        workspace.xH,
        workspace.dxH_dt,
        workspace.d2xH_dt,
        workspace.d3xH_dt,
        workspace.xBL,
        a,
        q,
        Vrr,
        ∂Vrr_∂t,
        workspace.Virr,
        workspace.∂Vrr_∂a,
        workspace.∂Virr_∂t,
        workspace.∂Virr_∂a;
        coordinates="harmonic",
    )

    copy_four_vector!(workspace.aSF25_H, workspace.aSF_H)
    copy_four_vector!(workspace.aSF25_BL, workspace.aSF_BL)
    iChimera.Inspiral.apply_bt_to_harmonic_rr_correction!(
        workspace.aSF25_H,
        workspace.aSF25_BL,
        workspace.aSF_corr_H,
        workspace.aSF_corr_BL,
        workspace.Acoord_corr_H,
        workspace.Acoord_BT,
        workspace.xH,
        workspace.dxH_dt,
        a,
        q;
        include25=true,
        include35=false,
    )

    copy_four_vector!(workspace.aSF35_H, workspace.aSF_H)
    copy_four_vector!(workspace.aSF35_BL, workspace.aSF_BL)
    iChimera.Inspiral.apply_bt_to_harmonic_rr_correction!(
        workspace.aSF35_H,
        workspace.aSF35_BL,
        workspace.aSF_corr_H,
        workspace.aSF_corr_BL,
        workspace.Acoord_corr_H,
        workspace.Acoord_BT,
        workspace.xH,
        workspace.dxH_dt,
        a,
        q;
        include25=true,
        include35=true,
    )

    Edot25, Ldot25, Qdot25 = local_fluxes_from_bl_four_accel(a, p, e, θmin, E, L, Q, C, rr, θθ, ϕϕ, dt_dτ, r_dot, θ_dot, ϕ_dot, workspace.aSF25_BL)
    Edot35, Ldot35, Qdot35 = local_fluxes_from_bl_four_accel(a, p, e, θmin, E, L, Q, C, rr, θθ, ϕϕ, dt_dτ, r_dot, θ_dot, ϕ_dot, workspace.aSF35_BL)

    return tt, Edot25, Ldot25, Qdot25, Edot35, Ldot35, Qdot35
end

function sample_local_chimera_bt_to_harmonic_fluxes(; a::Float64, p::Float64, e::Float64, θmin::Float64, sign_Lz::Int64,
                                                    q::Float64, psi0::Float64, chi0::Float64, phi0::Float64,
                                                    lmax_mass_fluxes::Int64, lmax_current_fluxes::Int64,
                                                    OnePN::Float64, TwoPN::Float64, TwoPointFivePN::Float64,
                                                    time_fraction_of_max_period::Float64, n_samples::Int,
                                                    reltol::Float64, abstol::Float64, maxiters::Int)
    n_samples >= 2 || error("n_samples must be at least 2.")

    t_geodesic, Ω, T_max = compute_geodesic_duration(a, p, e, θmin, sign_Lz, time_fraction_of_max_period)
    E, L, Q, C, _, p3, p4, zp, zm = iChimera.BLTimeGeodesics.compute_ODE_params(a, p, e, θmin, sign_Lz)

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

    workspace = allocate_rr_workspace()

    Edot25 = zeros(n_samples)
    Ldot25 = zeros(n_samples)
    Qdot25 = zeros(n_samples)
    Edot35 = zeros(n_samples)
    Ldot35 = zeros(n_samples)
    Qdot35 = zeros(n_samples)

    for (idx, t_sample) in enumerate(sample_times)
        if idx > 1
            step!(integrator, t_sample - integrator.t, true)
        end

        _, Edot25[idx], Ldot25[idx], Qdot25[idx], Edot35[idx], Ldot35[idx], Qdot35[idx] =
            compute_local_flux_pair!(
                workspace,
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
                lmax_mass_fluxes=lmax_mass_fluxes,
                lmax_current_fluxes=lmax_current_fluxes,
                OnePN=OnePN,
                TwoPN=TwoPN,
                TwoPointFivePN=TwoPointFivePN,
            )
    end

    return (
        times=sample_times,
        Edot25=Edot25,
        Ldot25=Ldot25,
        Qdot25=Qdot25,
        Edot35=Edot35,
        Ldot35=Ldot35,
        Qdot35=Qdot35,
        Ω=Ω,
        T_max=T_max,
        t_geodesic=t_geodesic,
        E=E,
        L=L,
        Q=Q,
        C=C,
    )
end

function diagnostic_text(a::Float64, p::Float64, e::Float64, q::Float64, θmin::Float64, result)
    return "a = $(a)\np = $(p)\ne = $(e)\nq = $(q)\nθmin = $(round(θmin; digits=4))\nT_max = $(round(result.T_max; digits=6)) M\nt_span = $(round(result.t_geodesic; digits=6)) M"
end

function plot_local_fluxes(result; a::Float64, p::Float64, e::Float64, q::Float64, θmin::Float64,
                           save_plots::Bool=true, results_dir::String=".", plot_tag::String="local_chimera_bt_to_harmonic_fluxes")
    mkpath(results_dir)

    labels = ["Chimera + 2.5PN corr.", "Chimera + 2.5PN + 3.5PN corr."]
    colors = [:royalblue, :tomato]
    linestyles = [:solid, :dash]
    linewidths = [3.0, 3.0]
    text = diagnostic_text(a, p, e, q, θmin, result)

    QuickPlots.plot11(
        [result.times, result.times],
        [result.Edot25, result.Edot35];
        colors=colors,
        labels=labels,
        linestyles=linestyles,
        linewidths=linewidths,
        xlabel=L"t / M",
        ylabel=L"\dot{E}",
        text=text,
        text_xloc=0.05,
        text_yloc=0.95,
        legend=true,
        position=:rb,
        save_plot=save_plots,
        fname=joinpath(results_dir, plot_tag * "_Edot.png"),
    )

    QuickPlots.plot11(
        [result.times, result.times],
        [result.Ldot25, result.Ldot35];
        colors=colors,
        labels=labels,
        linestyles=linestyles,
        linewidths=linewidths,
        xlabel=L"t / M",
        ylabel=L"\dot{L}_z",
        text=text,
        text_xloc=0.05,
        text_yloc=0.95,
        legend=true,
        position=:rb,
        save_plot=save_plots,
        fname=joinpath(results_dir, plot_tag * "_Ldot.png"),
    )

    QuickPlots.plot11(
        [result.times, result.times],
        [result.Qdot25, result.Qdot35];
        colors=colors,
        labels=labels,
        linestyles=linestyles,
        linewidths=linewidths,
        xlabel=L"t / M",
        ylabel=L"\dot{Q}",
        text=text,
        text_xloc=0.05,
        text_yloc=0.95,
        legend=true,
        position=:rb,
        save_plot=save_plots,
        fname=joinpath(results_dir, plot_tag * "_Qdot.png"),
    )

    return nothing
end

θmin = compute_theta_min_for_script(
    a,
    p,
    e,
    inclination_deg,
    inclination_type,
    sign_Lz;
    use_direct_theta_min=use_direct_theta_min,
    theta_min_direct=theta_min_direct,
)

flux_result = sample_local_chimera_bt_to_harmonic_fluxes(
    a=a,
    p=p,
    e=e,
    θmin=θmin,
    sign_Lz=sign_Lz,
    q=q,
    psi0=psi0,
    chi0=chi0,
    phi0=phi0,
    lmax_mass_fluxes=lmax_mass_fluxes,
    lmax_current_fluxes=lmax_current_fluxes,
    OnePN=OnePN,
    TwoPN=TwoPN,
    TwoPointFivePN=TwoPointFivePN,
    time_fraction_of_max_period=time_fraction_of_max_period,
    n_samples=n_samples,
    reltol=reltol,
    abstol=abstol,
    maxiters=maxiters,
)

plot_local_fluxes(
    flux_result;
    a=a,
    p=p,
    e=e,
    q=q,
    θmin=θmin,
    save_plots=save_plots,
    results_dir=results_dir,
    plot_tag=plot_tag,
)

println("Computed local hybrid Chimera BT->harmonic-correction fluxes on a fixed geodesic.")
println("θmin = $(θmin) rad")
println("Ω_r, Ω_θ, Ω_ϕ = $(flux_result.Ω)")
println("T_max = $(flux_result.T_max) M")
println("t_span = $(flux_result.t_geodesic) M")
if save_plots
    println("Saved plots to $(results_dir)")
end
