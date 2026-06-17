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

load_settings() = (600, 500, 1400, 600, 900, 900, 28, 30, 30, 24, 24, 25, true, true)

end

include(joinpath(@__DIR__, "..", "src", "QuickPlots.jl"))

const HRR = iChimera.HarmonicRRAcceleration

"""
User-editable parameters.

- Set `use_direct_theta_min = true` to bypass the inclination -> θmin mapping.
- `coordinates = "cartesian"` evaluates the same harmonic-PN local RR formula
  on the Cartesian coordinate state as a diagnostic comparison.
- `time_fraction_of_max_period` uses

      T_max = max(2π / Ω_r, 2π / Ω_θ)

  with BL-time frequencies, ignoring the formally infinite/degenerate special
  cases represented internally by very large placeholder values.
"""

# ----------------------------- EMRI / geodesic parameters -----------------------------
a = 0.5
p = 10.0
e = 0.2
q = 1.0e-5
coordinates = "harmonic"           # "harmonic" or "cartesian"

inclination_deg = 45.0
inclination_type = "iota"          # "iota" or "theta_inc"
sign_Lz = inclination_deg < 90.0 ? 1 : -1

use_direct_theta_min = false
theta_min_direct = π / 3           # used only if use_direct_theta_min = true

psi0 = 0.1
chi0 = 0.2
phi0 = 0.3

time_fraction_of_max_period = 10.0  # geodesic duration in units of max(T_r, T_θ)
n_samples = 400

reltol = 1.0e-13
abstol = 1.0e-13
maxiters = Int(1e8)

save_plots = true
results_dir = joinpath(@__DIR__, "..", "Results", "local_harmonic_pn_fluxes")
plot_tag = "a$(a)_p$(p)_e$(e)_q$(q)_$(coordinates)"
# --------------------------------------------------------------------------------------

function coordinate_module(coordinates::String)
    if coordinates == "harmonic"
        return iChimera.HarmonicCoords
    elseif coordinates == "cartesian"
        return iChimera.CartesianCoords
    end
    error("Unsupported coordinates=$(repr(coordinates)). Choose \"harmonic\" or \"cartesian\".")
end

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

function compute_local_flux_pair(integrator,
                                 a::Float64, p::Float64, e::Float64, θmin::Float64, q::Float64,
                                 E::Float64, L::Float64, Q::Float64, C::Float64,
                                 p3::Float64, p4::Float64, zp::Float64, zm::Float64;
                                 coordinates::String)
    tt, rr, θθ, ϕϕ, r_dot, θ_dot, ϕ_dot, _, _, _, dt_dτ, _, _ =
        iChimera.Inspiral.compute_geodesic_arrays(integrator, a, E, L, p, e, θmin, p3, p4, zp, zm)

    xBL = @SVector [rr, θθ, ϕϕ]
    vBL = @SVector [r_dot, θ_dot, ϕ_dot]

    coord_mod = coordinate_module(coordinates)
    xH = coord_mod.xBLtoH(xBL, a)
    vH = coord_mod.vBLtoH(xH, vBL, a)
    gH = coord_mod.g_μν_H(xH, a)

    Acoord25 = HRR.rr_harmonic_accel(xH, vH, q; include25=true, include35=false)
    Acoord35 = HRR.rr_harmonic_accel(xH, vH, q; include25=true, include35=true)

    a4H25 = HRR.lift_coord_accel_to_four(Acoord25, vH, gH)
    a4H35 = HRR.lift_coord_accel_to_four(Acoord35, vH, gH)

    a4BL25 = HRR.four_accel_coord_to_BL(a4H25, xH, a; coord_system=coordinates)
    a4BL35 = HRR.four_accel_coord_to_BL(a4H35, xH, a; coord_system=coordinates)

    Edot25, Ldot25, Qdot25 = local_fluxes_from_bl_four_accel(a, p, e, θmin, E, L, Q, C, rr, θθ, ϕϕ, dt_dτ, r_dot, θ_dot, ϕ_dot, a4BL25)
    Edot35, Ldot35, Qdot35 = local_fluxes_from_bl_four_accel(a, p, e, θmin, E, L, Q, C, rr, θθ, ϕϕ, dt_dτ, r_dot, θ_dot, ϕ_dot, a4BL35)

    return tt, Edot25, Ldot25, Qdot25, Edot35, Ldot35, Qdot35
end

function sample_local_harmonic_pn_fluxes(; a::Float64, p::Float64, e::Float64, θmin::Float64, sign_Lz::Int64,
                                         q::Float64, psi0::Float64, chi0::Float64, phi0::Float64,
                                         coordinates::String,
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
            compute_local_flux_pair(integrator, a, p, e, θmin, q, E, L, Q, C, p3, p4, zp, zm; coordinates=coordinates)
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

function diagnostic_text(a::Float64, p::Float64, e::Float64, q::Float64, θmin::Float64, coordinates::String, result)
    return "coords = $(coordinates)\na = $(a)\np = $(p)\ne = $(e)\nq = $(q)\nθmin = $(round(θmin; digits=4))\nT_max = $(round(result.T_max; digits=6)) M\nt_span = $(round(result.t_geodesic; digits=6)) M"
end

function get_standard_form_exponent(x::Float64)::Int64
    return iszero(x) ? 0 : floor(Int, log10(abs(x)))
end

function flux_scale_exponent(flux_series...)
    max_abs_flux = maximum(maximum(abs.(series)) for series in flux_series)
    return get_standard_form_exponent(max_abs_flux)
end

function plot_scaled_flux_pair(times::Vector{Float64}, flux25::Vector{Float64}, flux35::Vector{Float64};
                               ylabel,
                               save_plot::Bool,
                               fname::String)
    labels = ["2.5PN", "3.5PN"]
    colors = [:royalblue, :tomato]
    linestyles = [:solid, :dash]
    linewidths = [3.0, 3.0]

    scale_factor = flux_scale_exponent(flux25, flux35)
    scale_factor_label = L"\times 10^{%$(scale_factor)}"

    QuickPlots.with_theme(QuickPlots.theme_latexfonts()) do
        fig = QuickPlots.Figure(size = (QuickPlots.fig_width_1, QuickPlots.fig_height_1))
        ax = QuickPlots.Axis(
            fig[1, 1],
            ylabel = ylabel,
            xlabel = L"t / M",
            xlabelsize = QuickPlots.xlabelsize,
            ylabelsize = QuickPlots.ylabelsize,
            xticklabelsize = QuickPlots.xticklabelsize,
            yticklabelsize = QuickPlots.yticklabelsize,
            xgridvisible = QuickPlots.xgridvisible,
            ygridvisible = QuickPlots.ygridvisible,
        )

        QuickPlots.lines!(
            ax,
            times,
            flux25 .* 10.0^(-scale_factor),
            color = colors[1],
            linestyle = linestyles[1],
            linewidth = linewidths[1],
            label = labels[1],
        )
        QuickPlots.lines!(
            ax,
            times,
            flux35 .* 10.0^(-scale_factor),
            color = colors[2],
            linestyle = linestyles[2],
            linewidth = linewidths[2],
            label = labels[2],
        )

        QuickPlots.axislegend(ax, position = :rb, labelsize = QuickPlots.xticklabelsize, framevisible = true)
        QuickPlots.Label(fig[1, 1, QuickPlots.Top()], halign = :left, fontsize = QuickPlots.xticklabelsize, scale_factor_label)
        QuickPlots.resize_to_layout!(fig)
        save_plot && QuickPlots.save(fname, fig)
        display(fig)
    end

    return nothing
end

function plot_local_fluxes(result; a::Float64, p::Float64, e::Float64, q::Float64, θmin::Float64,
                           coordinates::String, save_plots::Bool=true, results_dir::String=".", plot_tag::String="local_harmonic_pn_fluxes")
    mkpath(results_dir)

    plot_scaled_flux_pair(
        result.times,
        result.Edot25,
        result.Edot35;
        ylabel = L"\dot{E}",
        save_plot = save_plots,
        fname = joinpath(results_dir, plot_tag * "_Edot.png"),
    )

    # plot_scaled_flux_pair(
    #     result.times,
    #     result.Ldot25,
    #     result.Ldot35;
    #     ylabel = L"\dot{L}_z",
    #     save_plot = save_plots,
    #     fname = joinpath(results_dir, plot_tag * "_Ldot.png"),
    # )

    # plot_scaled_flux_pair(
    #     result.times,
    #     result.Qdot25,
    #     result.Qdot35;
    #     ylabel = L"\dot{Q}",
    #     save_plot = save_plots,
    #     fname = joinpath(results_dir, plot_tag * "_Qdot.png"),
    # )

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

flux_result = sample_local_harmonic_pn_fluxes(
    a=a,
    p=p,
    e=e,
    θmin=θmin,
    sign_Lz=sign_Lz,
    q=q,
    psi0=psi0,
    chi0=chi0,
    phi0=phi0,
    coordinates=coordinates,
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
    coordinates=coordinates,
    save_plots=save_plots,
    results_dir=results_dir,
    plot_tag=plot_tag,
)

println("Computed local harmonic-PN fluxes on a fixed geodesic.")
println("coordinates = $(coordinates)")
println("θmin = $(θmin) rad")
println("Ω_r, Ω_θ, Ω_ϕ = $(flux_result.Ω)")
println("T_max = $(flux_result.T_max) M")
println("t_span = $(flux_result.t_geodesic) M")
if save_plots
    println("Saved plots to $(results_dir)")
end
