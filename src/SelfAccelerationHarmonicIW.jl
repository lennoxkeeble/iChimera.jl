module SelfAccelerationHarmonicIW

using LinearAlgebra
using StaticArrays

using ..HarmonicCoords
using ..SelfAccelerationHarmonic

const SH = SelfAccelerationHarmonic

################################################################################
# Iyer-Will / Burke-Thorne reaction-gauge transform for the harmonic Chimera
# self-force.
#
# Implements
#
#     hRR_{μν} -> hRR_{μν} - L_ξ gK_{μν}
#
# in harmonic coordinates.  The resulting general hRR_{μν} is then inserted into
# the same general MiSaTaQuWa-style expression used by Chimera:
#
#     a^α_RR = -Γ^2 P^{αβ} [ G^RR_{μνβ} u^μ u^ν
#                            - h^RR_{βγ} Γ^γ_{μν}[gK] u^μ u^ν ].
#
# This code assumes the same dimensionless units as iChimera: M_BH = 1.
# The EMRI mass ratio q is μ/M_BH, so η m^2 = q at leading EMRI order.
#
# The default current and target gauge are both Burke-Thorne:
#
#     α_current = 4, β_current = 5,
#     α_target  = 4, β_target  = 5.
#
# With these defaults there is no gauge change and the public aRRα function below
# calls the original SelfAccelerationHarmonic.aRRα exactly.  This makes the
# no-gauge-change regression test immediate.
################################################################################

# ---------------------------------------------------------------------------
# Basic helpers.
# Keep local versions here so this diagnostic module does not depend on whether
# SelfAccelerationHarmonic exports norm_3d/dot3d.
# ---------------------------------------------------------------------------

dot3d(u, v) = u[1]*v[1] + u[2]*v[2] + u[3]*v[3]
norm2_3d(u) = dot3d(u, u)
norm_3d(u) = sqrt(norm2_3d(u))

_fd_step(xH; epsrel=1.0e-5) = epsrel * max(1.0, norm_3d(xH))

_metricK(xH, spin) = Matrix{Float64}(HarmonicCoords.g_μν_H(xH, spin))
_metricKinv(xH, spin) = Matrix{Float64}(HarmonicCoords.gμν_H(xH, spin))

function _dgK_dx(xH, spin; eps_fd=_fd_step(xH))
    dg = zeros(4, 4, 3)

    for k in 1:3
        xp = collect(xH)
        xm = collect(xH)

        xp[k] += eps_fd
        xm[k] -= eps_fd

        dg[:, :, k] .= (_metricK(xp, spin) - _metricK(xm, spin)) ./ (2.0 * eps_fd)
    end

    return dg
end

_dg(dg, μ, ν, ρ) = ρ == 1 ? 0.0 : dg[μ, ν, ρ - 1]

function _christoffelK(xH, spin; eps_fd=_fd_step(xH))
    ginv = _metricKinv(xH, spin)
    dg = _dgK_dx(xH, spin; eps_fd=eps_fd)

    Γ = zeros(4, 4, 4)

    for α in 1:4, μ in 1:4, ν in 1:4, β in 1:4
        Γ[α, μ, ν] += 0.5 * ginv[α, β] * (
            _dg(dg, ν, β, μ) +
            _dg(dg, μ, β, ν) -
            _dg(dg, μ, ν, β)
        )
    end

    return Γ
end

# ---------------------------------------------------------------------------
# Build the original Chimera hRR_{μν} from the already-computed PM potentials:
#
#     h_tt = 2 Vrr,
#     h_ti = -4 Virr_i,
#     h_ij = 2 δ_ij Vrr.
#
# In fluxdev, RRPotentials.compute_RR_potentials! already returns Vrr and
# ∂Vrr_∂t and fills Virr, ∂Vrr_∂a, ∂Virr_∂t, ∂Virr_∂a.  Therefore this module
# should not recompute multipoles or potentials.
# ---------------------------------------------------------------------------

function _base_hRR(Vrr, Virr)
    h = zeros(4, 4)

    h[1, 1] = 2.0 * Vrr

    for i in 1:3
        h[1, i + 1] = -4.0 * Virr[i]
        h[i + 1, 1] = -4.0 * Virr[i]
        h[i + 1, i + 1] = 2.0 * Vrr
    end

    return h
end

function _base_dhRR(∂Vrr_∂t, ∂Vrr_∂a, ∂Virr_∂t, ∂Virr_∂a)
    dh = zeros(4, 4, 4)

    # Time derivative: coordinate index 1 corresponds to t.
    dh[1, 1, 1] = 2.0 * ∂Vrr_∂t

    for i in 1:3
        ip = i + 1

        dh[1, ip, 1] = -4.0 * ∂Virr_∂t[i]
        dh[ip, 1, 1] = -4.0 * ∂Virr_∂t[i]
        dh[ip, ip, 1] = 2.0 * ∂Vrr_∂t
    end

    # Spatial derivatives: coordinate indices 2,3,4 correspond to harmonic x,y,z.
    for k in 1:3
        kp = k + 1

        dh[1, 1, kp] = 2.0 * ∂Vrr_∂a[k]

        for i in 1:3
            ip = i + 1

            dh[1, ip, kp] = -4.0 * ∂Virr_∂a[i, k]
            dh[ip, 1, kp] = -4.0 * ∂Virr_∂a[i, k]
            dh[ip, ip, kp] = 2.0 * ∂Vrr_∂a[k]
        end
    end

    return dh
end

# ---------------------------------------------------------------------------
# Iyer-Will 2.5PN reaction-gauge vector.
#
# Sign convention:
#
#     x_target^i = x_current^i + ξ^i.
#
# The implementation uses
#
#     ξ^i = gauge_sign * (8q/15)/R * [ Δβ s n^i + C V^i ],
#
# with
#
#     Δα = α_target - α_current,
#     Δβ = β_target - β_current,
#     C  = 2Δβ - 3Δα.
#
# The default gauge_sign = -1.0 follows the convention used in the discussion.
# If a weak-field diagnostic gives the right magnitude but the opposite sign,
# flip gauge_sign -> +1.0.
# ---------------------------------------------------------------------------

function _gauge_data(xH, vH, aH, jerkH, q;
                     alpha_current=4.0,
                     beta_current=5.0,
                     alpha_target=4.0,
                     beta_target=5.0,
                     gauge_sign=-1.0)

    R = norm_3d(xH)
    R == 0.0 && error("Iyer-Will gauge vector is singular at R=0.")

    n = xH ./ R

    s = dot3d(n, vH)
    n_dot_a = dot3d(n, aH)
    n_dot_j = dot3d(n, jerkH)

    Δα = alpha_target - alpha_current
    Δβ = beta_target - beta_current
    C = 2.0 * Δβ - 3.0 * Δα

    # iChimera units have M_BH = 1.
    # At leading EMRI order η m^2 ≃ q.
    κ = 8.0 * q / 15.0
    pref = gauge_sign * κ

    ξ = zeros(3)
    dtξ = zeros(3)
    dttξ = zeros(3)
    dxξ = zeros(3, 3)
    dtdxξ = zeros(3, 3)

    for i in 1:3
        ξ[i] = pref / R * (Δβ * s * n[i] + C * vH[i])
        dtξ[i] = pref / R * (Δβ * n_dot_a * n[i] + C * aH[i])
        dttξ[i] = pref / R * (Δβ * n_dot_j * n[i] + C * jerkH[i])

        for j in 1:3
            δij = i == j ? 1.0 : 0.0

            dxξ[i, j] =
                pref / R^2 * (
                    Δβ * (n[i] * vH[j] + s * δij - 3.0 * s * n[i] * n[j])
                    - C * vH[i] * n[j]
                )

            dtdxξ[i, j] =
                pref / R^2 * (
                    Δβ * (n[i] * aH[j] + n_dot_a * δij - 3.0 * n_dot_a * n[i] * n[j])
                    - C * aH[i] * n[j]
                )
        end
    end

    return ξ, dtξ, dxξ, dttξ, dtdxξ
end

_ξderiv(dtξ, dxξ, k, μ) = μ == 1 ? dtξ[k] : dxξ[k, μ - 1]
_dt_ξderiv(dttξ, dtdxξ, k, μ) = μ == 1 ? dttξ[k] : dtdxξ[k, μ - 1]

# ---------------------------------------------------------------------------
# Lie derivative L_ξ gK in harmonic coordinates.
#
# Since ξ^t = 0 and the Kerr background is stationary:
#
#     (L_ξ g)_{μν}
#       = ξ^k ∂_k g_{μν}
#         + g_{kν} ∂_μ ξ^k
#         + g_{μk} ∂_ν ξ^k.
# ---------------------------------------------------------------------------

function _lie_xi_gK(xH, vH, aH, jerkH, q, spin;
                    eps_fd=_fd_step(xH),
                    alpha_current=4.0,
                    beta_current=5.0,
                    alpha_target=4.0,
                    beta_target=5.0,
                    gauge_sign=-1.0)

    g = _metricK(xH, spin)
    dg = _dgK_dx(xH, spin; eps_fd=eps_fd)

    ξ, dtξ, dxξ, _, _ = _gauge_data(
        xH, vH, aH, jerkH, q;
        alpha_current=alpha_current,
        beta_current=beta_current,
        alpha_target=alpha_target,
        beta_target=beta_target,
        gauge_sign=gauge_sign
    )

    L = zeros(4, 4)

    for μ in 1:4, ν in 1:4
        for k in 1:3
            ks = k + 1

            L[μ, ν] += ξ[k] * dg[μ, ν, k]
            L[μ, ν] += g[ks, ν] * _ξderiv(dtξ, dxξ, k, μ)
            L[μ, ν] += g[μ, ks] * _ξderiv(dtξ, dxξ, k, ν)
        end
    end

    return L
end

function _dt_lie_xi_gK(xH, vH, aH, jerkH, q, spin;
                       eps_fd=_fd_step(xH),
                       alpha_current=4.0,
                       beta_current=5.0,
                       alpha_target=4.0,
                       beta_target=5.0,
                       gauge_sign=-1.0)

    g = _metricK(xH, spin)
    dg = _dgK_dx(xH, spin; eps_fd=eps_fd)

    _, dtξ, _, dttξ, dtdxξ = _gauge_data(
        xH, vH, aH, jerkH, q;
        alpha_current=alpha_current,
        beta_current=beta_current,
        alpha_target=alpha_target,
        beta_target=beta_target,
        gauge_sign=gauge_sign
    )

    dtL = zeros(4, 4)

    for μ in 1:4, ν in 1:4
        for k in 1:3
            ks = k + 1

            dtL[μ, ν] += dtξ[k] * dg[μ, ν, k]
            dtL[μ, ν] += g[ks, ν] * _dt_ξderiv(dttξ, dtdxξ, k, μ)
            dtL[μ, ν] += g[μ, ks] * _dt_ξderiv(dttξ, dtdxξ, k, ν)
        end
    end

    return dtL
end

function _d_lie_xi_gK(xH, vH, aH, jerkH, q, spin;
                      eps_fd=_fd_step(xH),
                      alpha_current=4.0,
                      beta_current=5.0,
                      alpha_target=4.0,
                      beta_target=5.0,
                      gauge_sign=-1.0)

    dL = zeros(4, 4, 4)

    dL[:, :, 1] .= _dt_lie_xi_gK(
        xH, vH, aH, jerkH, q, spin;
        eps_fd=eps_fd,
        alpha_current=alpha_current,
        beta_current=beta_current,
        alpha_target=alpha_target,
        beta_target=beta_target,
        gauge_sign=gauge_sign
    )

    for k in 1:3
        xp = collect(xH)
        xm = collect(xH)

        xp[k] += eps_fd
        xm[k] -= eps_fd

        Lp = _lie_xi_gK(
            xp, vH, aH, jerkH, q, spin;
            eps_fd=eps_fd,
            alpha_current=alpha_current,
            beta_current=beta_current,
            alpha_target=alpha_target,
            beta_target=beta_target,
            gauge_sign=gauge_sign
        )

        Lm = _lie_xi_gK(
            xm, vH, aH, jerkH, q, spin;
            eps_fd=eps_fd,
            alpha_current=alpha_current,
            beta_current=beta_current,
            alpha_target=alpha_target,
            beta_target=beta_target,
            gauge_sign=gauge_sign
        )

        dL[:, :, k + 1] .= (Lp - Lm) ./ (2.0 * eps_fd)
    end

    return dL
end

# ---------------------------------------------------------------------------
# Gauge-transformed hRR and derivatives:
#
#     hRR -> hRR - L_ξ gK.
# ---------------------------------------------------------------------------

function _hRR_gauged(xH, vH, aH, jerkH, q, spin, Vrr, Virr;
                     eps_fd=_fd_step(xH),
                     alpha_current=4.0,
                     beta_current=5.0,
                     alpha_target=4.0,
                     beta_target=5.0,
                     gauge_sign=-1.0)

    h0 = _base_hRR(Vrr, Virr)

    Lξg = _lie_xi_gK(
        xH, vH, aH, jerkH, q, spin;
        eps_fd=eps_fd,
        alpha_current=alpha_current,
        beta_current=beta_current,
        alpha_target=alpha_target,
        beta_target=beta_target,
        gauge_sign=gauge_sign
    )

    return h0 - Lξg
end

function _dhRR_gauged(xH, vH, aH, jerkH, q, spin,
                      ∂Vrr_∂t, ∂Vrr_∂a, ∂Virr_∂t, ∂Virr_∂a;
                      eps_fd=_fd_step(xH),
                      alpha_current=4.0,
                      beta_current=5.0,
                      alpha_target=4.0,
                      beta_target=5.0,
                      gauge_sign=-1.0)

    dh0 = _base_dhRR(∂Vrr_∂t, ∂Vrr_∂a, ∂Virr_∂t, ∂Virr_∂a)

    dLξg = _d_lie_xi_gK(
        xH, vH, aH, jerkH, q, spin;
        eps_fd=eps_fd,
        alpha_current=alpha_current,
        beta_current=beta_current,
        alpha_target=alpha_target,
        beta_target=beta_target,
        gauge_sign=gauge_sign
    )

    return dh0 - dLξg
end

# ---------------------------------------------------------------------------
# Direct general-h_{μν} self-force.
#
# This avoids the scalar/vector-only Appendix-A shortcut, because a generic
# Iyer-Will gauge transform produces a non-isotropic spatial tensor h_ij.
# ---------------------------------------------------------------------------

function _Aβ_gauged(xH, vH, aH, jerkH, q, spin,
                    Vrr, ∂Vrr_∂t, Virr, ∂Vrr_∂a, ∂Virr_∂t, ∂Virr_∂a;
                    eps_fd=_fd_step(xH),
                    alpha_current=4.0,
                    beta_current=5.0,
                    alpha_target=4.0,
                    beta_target=5.0,
                    gauge_sign=-1.0)

    h = _hRR_gauged(
        xH, vH, aH, jerkH, q, spin, Vrr, Virr;
        eps_fd=eps_fd,
        alpha_current=alpha_current,
        beta_current=beta_current,
        alpha_target=alpha_target,
        beta_target=beta_target,
        gauge_sign=gauge_sign
    )

    dh = _dhRR_gauged(
        xH, vH, aH, jerkH, q, spin,
        ∂Vrr_∂t, ∂Vrr_∂a, ∂Virr_∂t, ∂Virr_∂a;
        eps_fd=eps_fd,
        alpha_current=alpha_current,
        beta_current=beta_current,
        alpha_target=alpha_target,
        beta_target=beta_target,
        gauge_sign=gauge_sign
    )

    ΓK = _christoffelK(xH, spin; eps_fd=eps_fd)

    u = [1.0, vH[1], vH[2], vH[3]]
    Aβ = zeros(4)

    for β in 1:4
        for μ in 1:4, ν in 1:4
            Gμνβ = 0.5 * (dh[ν, β, μ] + dh[μ, β, ν] - dh[μ, ν, β])
            Aβ[β] += Gμνβ * u[μ] * u[ν]

            for γ in 1:4
                Aβ[β] += -h[β, γ] * ΓK[γ, μ, ν] * u[μ] * u[ν]
            end
        end
    end

    return Aβ
end

# ---------------------------------------------------------------------------
# Public drop-in replacement for SelfAccelerationHarmonic.aRRα.
#
# If alpha_current == alpha_target and beta_current == beta_target, this calls
# the original fluxdev function exactly.  This makes the no-gauge-change test
# trivial and avoids finite-difference noise when no gauge transform is desired.
# ---------------------------------------------------------------------------

function aRRα(aSF_H,
              aSF_BL,
              xH,
              v,
              vH,
              aH,
              jerkH,
              xBL,
              rH,
              spin,
              q,
              Vrr,
              ∂Vrr_∂t,
              Virr,
              ∂Vrr_∂a,
              ∂Virr_∂t,
              ∂Virr_∂a;
              eps_fd=_fd_step(xH),
              alpha_current=4.0,
              beta_current=5.0,
              alpha_target=4.0,
              beta_target=5.0,
              gauge_sign=-1.0)

    if alpha_current == alpha_target && beta_current == beta_target
        SH.aRRα(
            aSF_H,
            aSF_BL,
            xH,
            v,
            vH,
            xBL,
            rH,
            spin,
            Vrr,
            ∂Vrr_∂t,
            Virr,
            ∂Vrr_∂a,
            ∂Virr_∂t,
            ∂Virr_∂a
        )
        return nothing
    end

    Aβ = _Aβ_gauged(
        xH,
        vH,
        aH,
        jerkH,
        q,
        spin,
        Vrr,
        ∂Vrr_∂t,
        Virr,
        ∂Vrr_∂a,
        ∂Virr_∂t,
        ∂Virr_∂a;
        eps_fd=eps_fd,
        alpha_current=alpha_current,
        beta_current=beta_current,
        alpha_target=alpha_target,
        beta_target=beta_target,
        gauge_sign=gauge_sign
    )

    aSF_H[:] = -SH.Γ(vH, xH, spin)^2 * SH.Pαβ(vH, xH, spin) * Aβ

    # Same convention as SelfAccelerationHarmonic.aRRα:
    # map the perturbing harmonic spatial acceleration back to BL using zero
    # velocity so that the Hessian/geodesic coordinate terms are not double-counted.
    aSF_BL[1] = aSF_H[1]
    aSF_BL[2:4] = HarmonicCoords.aHtoBL(xH, zeros(3), aSF_H[2:4], spin)

    return nothing
end

end # module SelfAccelerationHarmonicIW
