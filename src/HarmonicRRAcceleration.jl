module HarmonicRRAcceleration

using StaticArrays

using ..HarmonicCoords
using ..CartesianCoords

dot3d(u::AbstractVector, v::AbstractVector) = u[1] * v[1] + u[2] * v[2] + u[3] * v[3]
norm2_3d(u::AbstractVector) = dot3d(u, u)
norm_3d(u::AbstractVector) = sqrt(norm2_3d(u))

function _check_spatial_vector(name::AbstractString, u::AbstractVector)
    length(u) == 3 || throw(DimensionMismatch("$(name) must have length 3, got $(length(u))."))
    return nothing
end

function _check_four_vector(name::AbstractString, u::AbstractVector)
    length(u) == 4 || throw(DimensionMismatch("$(name) must have length 4, got $(length(u))."))
    return nothing
end

function _check_metric(gH::AbstractMatrix)
    size(gH) == (4, 4) || throw(DimensionMismatch("gH must be a 4x4 covariant metric, got size $(size(gH))."))
    return nothing
end

function _normalize_coord_system(coord_system::Symbol)
    if coord_system === :harmonic || coord_system === :cartesian
        return coord_system
    end
    throw(ArgumentError("Unsupported coord_system=$(repr(coord_system)). Supported values are :harmonic and :cartesian."))
end

_normalize_coord_system(coord_system::AbstractString) = _normalize_coord_system(Symbol(coord_system))

function _coord_module(coord_system::Union{Symbol, AbstractString})
    coord_sym = _normalize_coord_system(coord_system)
    return coord_sym === :harmonic ? HarmonicCoords : CartesianCoords
end

"""
    rr_harmonic_accel!(Acoord_H, X, V, q; M=1.0, G=1.0, c=1.0, include25=true, include35=true)

Fill `Acoord_H` with the harmonic-coordinate local radiation-reaction relative
coordinate acceleration

    d²X_H^i / dt² = a_2.5PN,H^i + a_3.5PN,H^i

from arXiv:2601.06743, specialized to the local 2.5PN + 3.5PN dissipative
terms. This is a spatial coordinate acceleration in Kerr-harmonic Cartesian
coordinates, not a four-acceleration and not the older Burke-Thorne-potential
Chimera force.
"""
function rr_harmonic_accel!(Acoord_H::AbstractVector, X::AbstractVector, V::AbstractVector, q::Real;
                            M::Real=1.0, G::Real=1.0, c::Real=1.0,
                            include25::Bool=true, include35::Bool=true)
    _check_spatial_vector("Acoord_H", Acoord_H)
    _check_spatial_vector("X", X)
    _check_spatial_vector("V", V)

    if !include25 && !include35
        Acoord_H[1] = zero(promote_type(eltype(Acoord_H), Float64))
        Acoord_H[2] = zero(promote_type(eltype(Acoord_H), Float64))
        Acoord_H[3] = zero(promote_type(eltype(Acoord_H), Float64))
        return Acoord_H
    end

    R2 = norm2_3d(X)
    R = sqrt(R2)
    R > 0 || throw(DomainError(R, "rr_harmonic_accel! is singular at R <= 0."))

    ν = q / (1 + q)^2
    m = M * (1 + q)
    u = G * m / R
    V2 = norm2_3d(V)
    Rdot = dot3d(X, V) / R

    invR = inv(R)
    invR3 = invR / R2

    N1 = X[1] * invR
    N2 = X[2] * invR
    N3 = X[3] * invR

    Acoord_H[1] = zero(promote_type(eltype(Acoord_H), typeof(float(R))))
    Acoord_H[2] = zero(promote_type(eltype(Acoord_H), typeof(float(R))))
    Acoord_H[3] = zero(promote_type(eltype(Acoord_H), typeof(float(R))))

    if include25
        pref25 = G^2 * m^2 * ν / (c^5 * R^3)
        coeffN25 = (136 / 15) * u * Rdot + (24 / 5) * Rdot * V2
        coeffV25 = (24 / 5) * u + (8 / 5) * V2

        Acoord_H[1] += pref25 * (N1 * coeffN25 - V[1] * coeffV25)
        Acoord_H[2] += pref25 * (N2 * coeffN25 - V[2] * coeffV25)
        Acoord_H[3] += pref25 * (N3 * coeffN25 - V[3] * coeffV25)
    end

    if include35
        pref35 = G^2 * m^2 * ν / (c^7 * R^3)

        NH =
            u^2 * Rdot * (-3956 / 35 - (184 / 5) * ν) +
            u * Rdot^3 * (-294 / 5 - (376 / 5) * ν) +
            u * Rdot * V2 * (-692 / 35 + (724 / 15) * ν) -
            112 * Rdot^5 +
            Rdot^3 * V2 * (114 + 12 * ν) +
            Rdot * V2^2 * (-366 / 35 - 12 * ν)

        VH =
            u^2 * (1060 / 21 + (104 / 5) * ν) +
            u * Rdot^2 * (82 / 3 + (848 / 15) * ν) +
            u * V2 * (-164 / 21 - (148 / 5) * ν) +
            120 * Rdot^4 +
            Rdot^2 * V2 * (-678 / 5 - (12 / 5) * ν) +
            V2^2 * (626 / 35 + (12 / 5) * ν)

        Acoord_H[1] += pref35 * (N1 * NH + V[1] * VH)
        Acoord_H[2] += pref35 * (N2 * NH + V[2] * VH)
        Acoord_H[3] += pref35 * (N3 * NH + V[3] * VH)
    end

    return Acoord_H
end

function rr_harmonic_accel(X::AbstractVector, V::AbstractVector, q::Real;
                           M::Real=1.0, G::Real=1.0, c::Real=1.0,
                           include25::Bool=true, include35::Bool=true)
    _check_spatial_vector("X", X)
    _check_spatial_vector("V", V)

    T = promote_type(typeof(float(X[1])), typeof(float(V[1])), typeof(float(q)), typeof(float(M)), typeof(float(G)), typeof(float(c)))
    Acoord_H = MVector{3, T}(undef)
    rr_harmonic_accel!(Acoord_H, X, V, q; M=M, G=G, c=c, include25=include25, include35=include35)
    return SVector{3, T}(Acoord_H)
end

"""
    rr_bt_accel!(Acoord_BT, X, V, q; M=1.0, G=1.0, c=1.0, include25=true, include35=true)

Fill `Acoord_BT` with the Burke-Thorne-gauge local radiation-reaction relative
coordinate acceleration

    d²X^i / dt² = a_2.5PN,BT^i + a_3.5PN,BT^i

through the same 2.5PN + 3.5PN local order used for the harmonic-coordinate
replacement model. This is a spatial coordinate acceleration, not a
four-acceleration. In the hybrid Chimera BT->harmonic model, this object is
used only inside the coordinate-acceleration correction

    ΔAcoord_H^i = Acoord_H^i - Acoord_BT^i .
"""
function rr_bt_accel!(Acoord_BT::AbstractVector, X::AbstractVector, V::AbstractVector, q::Real;
                      M::Real=1.0, G::Real=1.0, c::Real=1.0,
                      include25::Bool=true, include35::Bool=true)
    _check_spatial_vector("Acoord_BT", Acoord_BT)
    _check_spatial_vector("X", X)
    _check_spatial_vector("V", V)

    if !include25 && !include35
        Acoord_BT[1] = zero(promote_type(eltype(Acoord_BT), Float64))
        Acoord_BT[2] = zero(promote_type(eltype(Acoord_BT), Float64))
        Acoord_BT[3] = zero(promote_type(eltype(Acoord_BT), Float64))
        return Acoord_BT
    end

    R2 = norm2_3d(X)
    R = sqrt(R2)
    R > 0 || throw(DomainError(R, "rr_bt_accel! is singular at R <= 0."))

    ν = q / (1 + q)^2
    m = M * (1 + q)
    u = G * m / R
    V2 = norm2_3d(V)
    Rdot = dot3d(X, V) / R

    invR = inv(R)
    invR3 = invR / R2

    N1 = X[1] * invR
    N2 = X[2] * invR
    N3 = X[3] * invR

    Acoord_BT[1] = zero(promote_type(eltype(Acoord_BT), typeof(float(R))))
    Acoord_BT[2] = zero(promote_type(eltype(Acoord_BT), typeof(float(R))))
    Acoord_BT[3] = zero(promote_type(eltype(Acoord_BT), typeof(float(R))))

    if include25
        pref25 = 8 * G^2 * m^2 * ν / (5 * c^5 * R^3)
        A5 = 18 * V2 + (2 / 3) * u - 25 * Rdot^2
        B5 = 6 * V2 - 2 * u - 15 * Rdot^2

        Acoord_BT[1] += pref25 * (N1 * (A5 * Rdot) - V[1] * B5)
        Acoord_BT[2] += pref25 * (N2 * (A5 * Rdot) - V[2] * B5)
        Acoord_BT[3] += pref25 * (N3 * (A5 * Rdot) - V[3] * B5)
    end

    if include35
        pref35 = 8 * G^2 * m^2 * ν / (5 * c^7 * R^3)

        A7 =
            (87 / 14 - 48 * ν) * V2^2 -
            (5379 / 28 + (136 / 3) * ν) * V2 * u +
            (25 / 2) * (1 + 5 * ν) * V2 * Rdot^2 +
            (1353 / 4 + 133 * ν) * Rdot^2 * u -
            (35 / 2) * (1 - ν) * Rdot^4 +
            (160 / 7 + (55 / 3) * ν) * u^2

        B7 =
            -(27 / 14) * V2^2 -
            (4861 / 84 + (58 / 3) * ν) * V2 * u +
            (3 / 2) * (13 - 37 * ν) * V2 * Rdot^2 +
            (2591 / 12 + 97 * ν) * Rdot^2 * u -
            (25 / 2) * (1 - 7 * ν) * Rdot^4 +
            (1 / 3) * (776 / 7 + 55 * ν) * u^2

        Acoord_BT[1] += pref35 * (N1 * (A7 * Rdot) - V[1] * B7)
        Acoord_BT[2] += pref35 * (N2 * (A7 * Rdot) - V[2] * B7)
        Acoord_BT[3] += pref35 * (N3 * (A7 * Rdot) - V[3] * B7)
    end

    return Acoord_BT
end

function rr_bt_accel(X::AbstractVector, V::AbstractVector, q::Real;
                     M::Real=1.0, G::Real=1.0, c::Real=1.0,
                     include25::Bool=true, include35::Bool=true)
    _check_spatial_vector("X", X)
    _check_spatial_vector("V", V)

    T = promote_type(typeof(float(X[1])), typeof(float(V[1])), typeof(float(q)), typeof(float(M)), typeof(float(G)), typeof(float(c)))
    Acoord_BT = MVector{3, T}(undef)
    rr_bt_accel!(Acoord_BT, X, V, q; M=M, G=G, c=c, include25=include25, include35=include35)
    return SVector{3, T}(Acoord_BT)
end

"""
    rr_bt_to_harmonic_correction!(Acoord_corr_H, Acoord_BT, X, V, q; M=1.0, G=1.0, c=1.0, include25=true, include35=true)

Fill `Acoord_corr_H` with the local PN coordinate-acceleration correction

    ΔAcoord_H^i = Acoord_H^i - Acoord_BT^i

between the harmonic-coordinate and Burke-Thorne-coordinate relative RR
accelerations through 2.5PN + 3.5PN order. `Acoord_BT` is scratch storage used
to avoid allocations in hot loops.

This correction is the ingredient used in the hybrid Chimera BT->harmonic
branch. It is a coordinate-acceleration correction, not by itself a full
worldline gauge transformation.
"""
function rr_bt_to_harmonic_correction!(Acoord_corr_H::AbstractVector, Acoord_BT::AbstractVector, X::AbstractVector, V::AbstractVector, q::Real;
                                       M::Real=1.0, G::Real=1.0, c::Real=1.0,
                                       include25::Bool=true, include35::Bool=true)
    _check_spatial_vector("Acoord_corr_H", Acoord_corr_H)
    _check_spatial_vector("Acoord_BT", Acoord_BT)

    rr_harmonic_accel!(Acoord_corr_H, X, V, q; M=M, G=G, c=c, include25=include25, include35=include35)
    rr_bt_accel!(Acoord_BT, X, V, q; M=M, G=G, c=c, include25=include25, include35=include35)

    @inbounds for i in 1:3
        Acoord_corr_H[i] -= Acoord_BT[i]
    end

    return Acoord_corr_H
end

function rr_bt_to_harmonic_correction(X::AbstractVector, V::AbstractVector, q::Real;
                                      M::Real=1.0, G::Real=1.0, c::Real=1.0,
                                      include25::Bool=true, include35::Bool=true)
    _check_spatial_vector("X", X)
    _check_spatial_vector("V", V)

    T = promote_type(typeof(float(X[1])), typeof(float(V[1])), typeof(float(q)), typeof(float(M)), typeof(float(G)), typeof(float(c)))
    Acoord_corr_H = MVector{3, T}(undef)
    Acoord_BT = MVector{3, T}(undef)
    rr_bt_to_harmonic_correction!(Acoord_corr_H, Acoord_BT, X, V, q; M=M, G=G, c=c, include25=include25, include35=include35)
    return SVector{3, T}(Acoord_corr_H)
end

"""
    rr_bt_to_harmonic_correction_four_accel_BL!(a4_corr_BL, a4_corr_H, Acoord_corr_H, Acoord_BT, X, V, spin, q; ...)

Fill `a4_corr_BL` with the approach-2 correction

    δa_BL = lift_to_BL(a_PN,H - a_PN,BT)

using the harmonic-coordinate position `X` and coordinate velocity `V`. The
scratch arrays `a4_corr_H`, `Acoord_corr_H`, and `Acoord_BT` are supplied by the
caller to avoid allocations in local-flux loops. This function computes only
the additive correction; it does not recompute the Chimera contribution.
"""
function rr_bt_to_harmonic_correction_four_accel_BL!(a4_corr_BL::AbstractVector,
                                                     a4_corr_H::AbstractVector,
                                                     Acoord_corr_H::AbstractVector,
                                                     Acoord_BT::AbstractVector,
                                                     X::AbstractVector,
                                                     V::AbstractVector,
                                                     spin::Real,
                                                     q::Real;
                                                     M::Real=1.0,
                                                     G::Real=1.0,
                                                     c::Real=1.0,
                                                     include25::Bool=true,
                                                     include35::Bool=true)
    _check_four_vector("a4_corr_BL", a4_corr_BL)
    _check_four_vector("a4_corr_H", a4_corr_H)
    _check_spatial_vector("Acoord_corr_H", Acoord_corr_H)
    _check_spatial_vector("Acoord_BT", Acoord_BT)
    _check_spatial_vector("X", X)
    _check_spatial_vector("V", V)

    rr_bt_to_harmonic_correction!(Acoord_corr_H, Acoord_BT, X, V, q; M=M, G=G, c=c, include25=include25, include35=include35)
    gH = HarmonicCoords.g_μν_H(X, spin)
    lift_coord_accel_to_four!(a4_corr_H, Acoord_corr_H, V, gH)
    four_accel_H_to_BL!(a4_corr_BL, a4_corr_H, X, spin)
    return a4_corr_BL
end

function rr_bt_to_harmonic_correction_four_accel_BL(X::AbstractVector,
                                                    V::AbstractVector,
                                                    spin::Real,
                                                    q::Real;
                                                    M::Real=1.0,
                                                    G::Real=1.0,
                                                    c::Real=1.0,
                                                    include25::Bool=true,
                                                    include35::Bool=true)
    _check_spatial_vector("X", X)
    _check_spatial_vector("V", V)

    T = promote_type(typeof(float(X[1])), typeof(float(V[1])), typeof(float(spin)), typeof(float(q)), typeof(float(M)), typeof(float(G)), typeof(float(c)))
    a4_corr_BL = MVector{4, T}(undef)
    a4_corr_H = MVector{4, T}(undef)
    Acoord_corr_H = MVector{3, T}(undef)
    Acoord_BT = MVector{3, T}(undef)

    rr_bt_to_harmonic_correction_four_accel_BL!(
        a4_corr_BL,
        a4_corr_H,
        Acoord_corr_H,
        Acoord_BT,
        X,
        V,
        spin,
        q;
        M=M,
        G=G,
        c=c,
        include25=include25,
        include35=include35,
    )
    return SVector{4, T}(a4_corr_BL)
end

"""
    add_bt_to_harmonic_correction_four_accel_BL!(a4_fixed_BL, a4_chimera_BL, ...)

Fill `a4_fixed_BL` with `a4_chimera_BL + δa_BL`, where `δa_BL` is the
approach-2 BT-to-harmonic correction. This is intended for post-processing
local Chimera fluxes without recomputing the Chimera self-force.
"""
function add_bt_to_harmonic_correction_four_accel_BL!(a4_fixed_BL::AbstractVector,
                                                      a4_chimera_BL::AbstractVector,
                                                      a4_corr_BL::AbstractVector,
                                                      a4_corr_H::AbstractVector,
                                                      Acoord_corr_H::AbstractVector,
                                                      Acoord_BT::AbstractVector,
                                                      X::AbstractVector,
                                                      V::AbstractVector,
                                                      spin::Real,
                                                      q::Real;
                                                      M::Real=1.0,
                                                      G::Real=1.0,
                                                      c::Real=1.0,
                                                      include25::Bool=true,
                                                      include35::Bool=true)
    _check_four_vector("a4_fixed_BL", a4_fixed_BL)
    _check_four_vector("a4_chimera_BL", a4_chimera_BL)

    rr_bt_to_harmonic_correction_four_accel_BL!(
        a4_corr_BL,
        a4_corr_H,
        Acoord_corr_H,
        Acoord_BT,
        X,
        V,
        spin,
        q;
        M=M,
        G=G,
        c=c,
        include25=include25,
        include35=include35,
    )

    @inbounds for i in 1:4
        a4_fixed_BL[i] = a4_chimera_BL[i] + a4_corr_BL[i]
    end

    return a4_fixed_BL
end

"""
    gamma_from_metric_velocity(gH, V_H)

Return `u^t = Γ` from the covariant harmonic-coordinate metric `gH = g_{μν}`
and the coordinate velocity `V_H = dX_H/dt`. This assumes the repo convention
with signature `(-,+,+,+)`.
"""
function gamma_from_metric_velocity(gH::AbstractMatrix, V_H::AbstractVector)
    _check_metric(gH)
    _check_spatial_vector("V_H", V_H)

    timelike_norm =
        gH[1, 1] +
        2 * (gH[1, 2] * V_H[1] + gH[1, 3] * V_H[2] + gH[1, 4] * V_H[3]) +
        gH[2, 2] * V_H[1]^2 +
        gH[3, 3] * V_H[2]^2 +
        gH[4, 4] * V_H[3]^2 +
        2 * gH[2, 3] * V_H[1] * V_H[2] +
        2 * gH[2, 4] * V_H[1] * V_H[3] +
        2 * gH[3, 4] * V_H[2] * V_H[3]

    denom = -timelike_norm
    denom > 0 || throw(DomainError(denom, "gamma_from_metric_velocity expected a timelike coordinate velocity with -g(u,u) > 0."))
    return inv(sqrt(denom))
end

"""
    coord_accel_from_four!(Acoord_H, a4_H, V_H, Γ)

Recover the harmonic coordinate acceleration from a contravariant
four-acceleration and coordinate velocity via

    Acoord_H^i = Γ^{-2} (a^i - V^i a^t).
"""
function coord_accel_from_four!(Acoord_H::AbstractVector, a4_H::AbstractVector, V_H::AbstractVector, Γ::Real)
    _check_spatial_vector("Acoord_H", Acoord_H)
    _check_four_vector("a4_H", a4_H)
    _check_spatial_vector("V_H", V_H)

    invΓ2 = inv(Γ^2)
    Acoord_H[1] = (a4_H[2] - V_H[1] * a4_H[1]) * invΓ2
    Acoord_H[2] = (a4_H[3] - V_H[2] * a4_H[1]) * invΓ2
    Acoord_H[3] = (a4_H[4] - V_H[3] * a4_H[1]) * invΓ2
    return Acoord_H
end

function coord_accel_from_four(a4_H::AbstractVector, V_H::AbstractVector, Γ::Real)
    T = promote_type(typeof(float(a4_H[1])), typeof(float(V_H[1])), typeof(float(Γ)))
    Acoord_H = MVector{3, T}(undef)
    coord_accel_from_four!(Acoord_H, a4_H, V_H, Γ)
    return SVector{3, T}(Acoord_H)
end

"""
    lift_coord_accel_to_four!(a4_H, Acoord_H, V_H, gH)

Lift the harmonic coordinate acceleration `Acoord_H = d²X_H/dt²` to the
contravariant harmonic-coordinate four-acceleration `a4_H` using the covariant
metric `gH = g_{μν}` and the orthogonality relation `u_μ a^μ = 0`.
"""
function lift_coord_accel_to_four!(a4_H::AbstractVector, Acoord_H::AbstractVector, V_H::AbstractVector, gH::AbstractMatrix)
    _check_four_vector("a4_H", a4_H)
    _check_spatial_vector("Acoord_H", Acoord_H)
    _check_spatial_vector("V_H", V_H)
    _check_metric(gH)

    Γ = gamma_from_metric_velocity(gH, V_H)

    u1 = Γ
    u2 = Γ * V_H[1]
    u3 = Γ * V_H[2]
    u4 = Γ * V_H[3]

    u_cov_1 = gH[2, 1] * u1 + gH[2, 2] * u2 + gH[2, 3] * u3 + gH[2, 4] * u4
    u_cov_2 = gH[3, 1] * u1 + gH[3, 2] * u2 + gH[3, 3] * u3 + gH[3, 4] * u4
    u_cov_3 = gH[4, 1] * u1 + gH[4, 2] * u2 + gH[4, 3] * u3 + gH[4, 4] * u4

    a_t = Γ^3 * (u_cov_1 * Acoord_H[1] + u_cov_2 * Acoord_H[2] + u_cov_3 * Acoord_H[3])

    a4_H[1] = a_t
    a4_H[2] = Γ^2 * Acoord_H[1] + V_H[1] * a_t
    a4_H[3] = Γ^2 * Acoord_H[2] + V_H[2] * a_t
    a4_H[4] = Γ^2 * Acoord_H[3] + V_H[3] * a_t

    return a4_H
end

function lift_coord_accel_to_four(Acoord_H::AbstractVector, V_H::AbstractVector, gH::AbstractMatrix)
    T = promote_type(typeof(float(Acoord_H[1])), typeof(float(V_H[1])), typeof(float(gH[1, 1])))
    a4_H = MVector{4, T}(undef)
    lift_coord_accel_to_four!(a4_H, Acoord_H, V_H, gH)
    return SVector{4, T}(a4_H)
end

"""
    four_accel_H_to_BL!(a4_BL, a4_H, xH, spin)

Transform a contravariant harmonic-coordinate four-acceleration into BL
components. The spatial part is transformed as a four-vector with the harmonic
to BL Jacobian only. No Hessian terms appear here because the input is already a
four-vector, not a coordinate acceleration.
"""
function four_accel_H_to_BL!(a4_BL::AbstractVector, a4_H::AbstractVector, xH::AbstractVector, spin::Real)
    _check_four_vector("a4_BL", a4_BL)
    _check_four_vector("a4_H", a4_H)
    _check_spatial_vector("xH", xH)

    jBLH = HarmonicCoords.jBLH(xH, spin)

    a4_BL[1] = a4_H[1]
    a4_BL[2] = jBLH[1, 1] * a4_H[2] + jBLH[1, 2] * a4_H[3] + jBLH[1, 3] * a4_H[4]
    a4_BL[3] = jBLH[2, 1] * a4_H[2] + jBLH[2, 2] * a4_H[3] + jBLH[2, 3] * a4_H[4]
    a4_BL[4] = jBLH[3, 1] * a4_H[2] + jBLH[3, 2] * a4_H[3] + jBLH[3, 3] * a4_H[4]

    return a4_BL
end

function four_accel_H_to_BL(a4_H::AbstractVector, xH::AbstractVector, spin::Real)
    T = promote_type(typeof(float(a4_H[1])), typeof(float(xH[1])), typeof(float(spin)))
    a4_BL = MVector{4, T}(undef)
    four_accel_H_to_BL!(a4_BL, a4_H, xH, spin)
    return SVector{4, T}(a4_BL)
end

"""
    four_accel_coord_to_BL!(a4_BL, a4_coord, xcoord, spin; coord_system=:harmonic)

Transform a contravariant four-acceleration from either Kerr-harmonic or the
repo's Cartesian coordinate chart into BL components. The spatial part is
transformed as a four-vector with the appropriate Jacobian only.
"""
function four_accel_coord_to_BL!(a4_BL::AbstractVector, a4_coord::AbstractVector, xcoord::AbstractVector, spin::Real;
                                 coord_system::Union{Symbol, AbstractString}=:harmonic)
    _check_four_vector("a4_BL", a4_BL)
    _check_four_vector("a4_coord", a4_coord)
    _check_spatial_vector("xcoord", xcoord)

    coord_mod = _coord_module(coord_system)
    jBLcoord = coord_mod.jBLH(xcoord, spin)

    a4_BL[1] = a4_coord[1]
    a4_BL[2] = jBLcoord[1, 1] * a4_coord[2] + jBLcoord[1, 2] * a4_coord[3] + jBLcoord[1, 3] * a4_coord[4]
    a4_BL[3] = jBLcoord[2, 1] * a4_coord[2] + jBLcoord[2, 2] * a4_coord[3] + jBLcoord[2, 3] * a4_coord[4]
    a4_BL[4] = jBLcoord[3, 1] * a4_coord[2] + jBLcoord[3, 2] * a4_coord[3] + jBLcoord[3, 3] * a4_coord[4]

    return a4_BL
end

function four_accel_coord_to_BL(a4_coord::AbstractVector, xcoord::AbstractVector, spin::Real;
                                coord_system::Union{Symbol, AbstractString}=:harmonic)
    T = promote_type(typeof(float(a4_coord[1])), typeof(float(xcoord[1])), typeof(float(spin)))
    a4_BL = MVector{4, T}(undef)
    four_accel_coord_to_BL!(a4_BL, a4_coord, xcoord, spin; coord_system=coord_system)
    return SVector{4, T}(a4_BL)
end

end
