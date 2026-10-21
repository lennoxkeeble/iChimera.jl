module RRPotentials
using StaticArrays

δ(x::Int, y::Int)::Int = x == y ? 1 : 0

const levi_civita_table = Dict(
    (1, 2, 3) => 1,
    (2, 3, 1) => 1,
    (3, 1, 2) => 1,
    (3, 2, 1) => -1,
    (2, 1, 3) => -1,
    (1, 3, 2) => -1
)

ε(i::Int, j::Int, k::Int)::Int = get(levi_civita_table, (i, j, k), 0)

function compute_RR_potentials!(
    Virr::MVector{3, Float64},
    ∂Vrr_∂a::MVector{3, Float64},
    ∂Virr_∂t::MVector{3, Float64},
    ∂Virr_∂a::MMatrix{3, 3, Float64, 9},
    x::MVector{3, Float64},
    Mij5::MMatrix{3, 3, Float64, 9},
    Mij6::MMatrix{3, 3, Float64, 9},
    Mij7::MMatrix{3, 3, Float64, 9},
    Mij8::MMatrix{3, 3, Float64, 9},
    Mijk7::MArray{Tuple{3, 3, 3}, Float64, 3, 27},
    Mijk8::MArray{Tuple{3, 3, 3}, Float64, 3, 27},
    Skl5::MMatrix{3, 3, Float64, 9},
    Skl6::MMatrix{3, 3, Float64, 9};
    ΔMij5::Union{Nothing, MMatrix{3, 3, Float64, 9}} = nothing,
    ΔMij6::Union{Nothing, MMatrix{3, 3, Float64, 9}} = nothing,
)
    Vrr_value = Vrr(x, Mij5, Mij7, Mijk7)
    ∂Vrr_∂t_value = ∂Vrr_∂t(x, Mij6, Mij8, Mijk8)

    for i in 1:3
        Virr[i] = RRPotentials.Virr(i, x, Mij6, Skl5)
        ∂Vrr_∂a[i] = ∂Vrr_∂xm(i, x, Mij5, Mij7, Mijk7)
        ∂Virr_∂t[i] = RRPotentials.∂Virr_∂t(i, x, Mij7, Skl6)
        for j in 1:3
            ∂Virr_∂a[i, j] = ∂Virr_∂xm(i, j, x, Mij6, Skl5)
        end
    end

    # optional 1.5pn tail (blanchet 1997, eq. 4.33a): a burke-thorne term with M^(5), M^(6) replaced by the tail integrals
    # ΔMij5, ΔMij6 (see tail_weights); V_i is unchanged at this order
    (ΔMij5 === nothing) == (ΔMij6 === nothing) || throw(ArgumentError("ΔMij5 and ΔMij6 must be given together"))
    if ΔMij5 !== nothing
        zero_ij = zero(MMatrix{3, 3, Float64, 9})
        zero_ijk = zero(MArray{Tuple{3, 3, 3}, Float64, 3, 27})
        Vrr_value += Vrr(x, ΔMij5, zero_ij, zero_ijk)
        ∂Vrr_∂t_value += ∂Vrr_∂t(x, ΔMij6, zero_ij, zero_ijk)
        for i in 1:3
            ∂Vrr_∂a[i] += ∂Vrr_∂xm(i, x, ΔMij5, zero_ij, zero_ijk)
        end
    end

    return Vrr_value, ∂Vrr_∂t_value
end

# quadrature weights for the 1.5pn tail integral of blanchet 1997, eq. 4.33a,
#   ΔM(t_n) = 4M ∫_0^∞ dλ ln(λ/2) f(t_n - λ) ≈ Σ_{j=0}^{n_past} W[j+1] f(t_n - j h),
# with f = M^(7) for ΔMij5 and f = M^(8) for ΔMij6, sampled uniformly with step h (time scale of the logarithm: 1 in units of M).
# the past is damped by exp(-ε λ) for each ε in ε_values and the results are extrapolated polynomially to ε = 0, which defines
# the integral for a source that oscillates indefinitely (e.g. a geodesic continued into the past). n_past must be even and
# exp(-minimum(ε_values) n_past h) negligible. quadrature: ln(λ/2) times a quadratic integrated exactly on [0, 2h], composite
# simpson rule on [2h, n_past h].
function tail_weights(h::Float64, n_past::Int64, ε_values::AbstractVector{Float64}; M::Float64 = 1.0)::Vector{Float64}
    (iseven(n_past) && n_past >= 2) || throw(ArgumentError("n_past must be even and at least 2"))
    λ = h .* (0:n_past)
    base = zeros(n_past + 1)
    # ∫_0^2 ℓ_k(u) du and ∫_0^2 ln(u) ℓ_k(u) du for the quadratic lagrange basis on u = 0, 1, 2
    simpson = (1 / 3, 4 / 3, 1 / 3)
    log_moments = (log(2) / 3 - 17 / 18, 4 * log(2) / 3 - 10 / 9, log(2) / 3 + 1 / 18)
    for k in 1:3
        base[k] += h * (log(h / 2) * simpson[k] + log_moments[k])
    end
    for j in 2:n_past
        w = (j == 2 || j == n_past) ? 1 / 3 : (isodd(j) ? 4 / 3 : 2 / 3)
        base[j + 1] += h * w * log(λ[j + 1] / 2)
    end
    W = zeros(n_past + 1)
    for (k, εk) in enumerate(ε_values)
        c = prod((-εl / (εk - εl) for (l, εl) in enumerate(ε_values) if l != k); init = 1.0)
        W .+= c .* exp.(-εk .* λ) .* base
    end
    return 4.0 * M .* W
end

function Vrr(
    x::MVector{3, Float64},
    Mij5::MMatrix{3, 3, Float64, 9},
    Mij7::MMatrix{3, 3, Float64, 9},
    Mijk7::MArray{Tuple{3, 3, 3}, Float64, 3, 27},
)::Float64
    value = 0.0
    r2 = sum(abs2, x)
    for i in 1:3, j in 1:3
        value -= x[i] * x[j] * Mij5[i, j] / 5.0
        value -= r2 * x[i] * x[j] * Mij7[i, j] / 70.0
        for k in 1:3
            value += x[i] * x[j] * x[k] * Mijk7[i, j, k] / 189.0
        end
    end
    return value
end

function ∂Vrr_∂t(
    x::MVector{3, Float64},
    Mij6::MMatrix{3, 3, Float64, 9},
    Mij8::MMatrix{3, 3, Float64, 9},
    Mijk8::MArray{Tuple{3, 3, 3}, Float64, 3, 27},
)::Float64
    value = 0.0
    r2 = sum(abs2, x)
    for i in 1:3, j in 1:3
        value -= x[i] * x[j] * Mij6[i, j] / 5.0
        value -= r2 * x[i] * x[j] * Mij8[i, j] / 70.0
        for k in 1:3
            value += x[i] * x[j] * x[k] * Mijk8[i, j, k] / 189.0
        end
    end
    return value
end

function ∂Vrr_∂xm(
    m::Int64,
    x::MVector{3, Float64},
    Mij5::MMatrix{3, 3, Float64, 9},
    Mij7::MMatrix{3, 3, Float64, 9},
    Mijk7::MArray{Tuple{3, 3, 3}, Float64, 3, 27},
)::Float64
    value = 0.0
    r2 = sum(abs2, x)
    for i in 1:3, j in 1:3
        value -= Mij5[i, j] * (x[j] * δ(i, m) + x[i] * δ(j, m)) / 5.0
        value -= Mij7[i, j] * (
            2.0 * x[m] * x[i] * x[j] +
            r2 * (x[j] * δ(i, m) + x[i] * δ(j, m))
        ) / 70.0
        for k in 1:3
            value += Mijk7[i, j, k] * (
                x[j] * x[k] * δ(i, m) +
                x[i] * x[k] * δ(j, m) +
                x[i] * x[j] * δ(k, m)
            ) / 189.0
        end
    end
    return value
end

function Virr(
    i::Int64,
    x::MVector{3, Float64},
    Mij6::MMatrix{3, 3, Float64, 9},
    Skl5::MMatrix{3, 3, Float64, 9},
)::Float64
    value = 0.0
    r2 = sum(abs2, x)
    for j in 1:3, k in 1:3
        value += (
            x[i] * x[j] * x[k] -
            r2 * (x[k] * δ(i, j) + x[j] * δ(i, k) + x[i] * δ(j, k)) / 5.0
        ) * Mij6[j, k] / 21.0
        for l in 1:3
            value -= 4.0 * x[j] * x[l] * ε(i, j, k) * Skl5[k, l] / 45.0
        end
    end
    return value
end

function ∂Virr_∂t(
    i::Int64,
    x::MVector{3, Float64},
    Mij7::MMatrix{3, 3, Float64, 9},
    Skl6::MMatrix{3, 3, Float64, 9},
)::Float64
    value = 0.0
    r2 = sum(abs2, x)
    for j in 1:3, k in 1:3
        value += (
            x[i] * x[j] * x[k] -
            r2 * (x[k] * δ(i, j) + x[j] * δ(i, k) + x[i] * δ(j, k)) / 5.0
        ) * Mij7[j, k] / 21.0
        for l in 1:3
            value -= 4.0 * x[j] * x[l] * ε(i, j, k) * Skl6[k, l] / 45.0
        end
    end
    return value
end

function ∂Virr_∂xm(
    i::Int64,
    m::Int64,
    x::MVector{3, Float64},
    Mij6::MMatrix{3, 3, Float64, 9},
    Skl5::MMatrix{3, 3, Float64, 9},
)::Float64
    value = 0.0
    r2 = sum(abs2, x)
    for j in 1:3, k in 1:3
        value += Mij6[j, k] * (
            x[j] * x[k] * δ(i, m) +
            x[i] * x[k] * δ(j, m) +
            x[i] * x[j] * δ(k, m) -
            (
                2.0 * x[m] * (x[k] * δ(i, j) + x[j] * δ(i, k) + x[i] * δ(j, k)) +
                r2 * (δ(i, m) * δ(j, k) + δ(i, k) * δ(j, m) + δ(i, j) * δ(k, m))
            ) / 5.0
        ) / 21.0
        for l in 1:3
            value -= 4.0 * Skl5[k, l] * (
                x[l] * δ(j, m) + x[j] * δ(l, m)
            ) * ε(i, j, k) / 45.0
        end
    end
    return value
end

end
