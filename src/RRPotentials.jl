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
    Skl6::MMatrix{3, 3, Float64, 9},
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

    return Vrr_value, ∂Vrr_∂t_value
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
