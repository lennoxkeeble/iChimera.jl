#= 

    In this module we write the analytic expressions for computing the self-force from the high order derivatives of the multipole moments and the position, velocity and acceleration in harmonic coordinates. 
    See Eqs. 54-56, 57, 61-63, and A1-A14 in arXiv:1109.0572v2. Note that our implemenation of Eqs. A12 - A14 differ slightly by factors of 2 in places and symmetrization operators.

=#

module SelfAccelerationHarmonic2
using LinearAlgebra
using Combinatorics
using ..HarmonicCoords
using StaticArrays
using ..RRPotentials
import ..Kerr.KerrMetric: g_μν, Γαμν

"""
# Common Arguments in this module
- `xBL::AbstractVector{Float64}`: Boyer-Lindquist coordinates, xBL = [r, θ, ϕ].
- `xH::AbstractVector{Float64}`: Harmonic coordinates, xH = [x, y, z].
- `vH::AbstractArray`: velocity in Harmonic coordinates.
- `v::Float64`: velocity v = sqrt(vx^2 + vy^2 + vz^2).
- `aH::AbstractArray`: acceleration in Harmonic coordinates.
- `rH::Float64`: rH = sqrt(xH^2 + yH^2 + zH^2).
- `jBLH::AbstractArray`: Jacobian of the transformation from BL to Harmonic coordinates.
- `HessBLH::AbstractArray`: Hessian of the transformation from BL to Harmonic coordinates.
- `Mij5::AbstractArray`: fifth derivative of the mass quadrupole (Eq. 48).
- `Mij6::AbstractArray`: sixth derivative of the mass quadrupole (Eq. 48).
- `Mij7::AbstractArray`: seventh derivative of the mass quadrupole (Eq. 48).
- `Mij8::AbstractArray`: eighth derivative of the mass quadrupole (Eq. 48).
- `Mijk7::AbstractArray`: seventh derivative of the mass quadrupole (Eq. 48).
- `Mijk8::AbstractArray`: eighth derivative of the mass quadrupole (Eq. 48).
- `Sij5::AbstractArray`: fifth derivative of the current quadrupole (Eq. 49).
- `Sij6::AbstractArray`: sixth derivative of the current quadrupole (Eq. 49).
- `∂Vrr_∂t::Float64`: time derivative of the radiation reaction potential (Eq. 44).
- `∂Vrr_∂a::AbstractVector{Float64}`: radiation reaction potential derivative with respect to the harmonic spatial coordinates.
- `∂Virr_∂t::AbstractVector{Float64}`: time derivative of the spatial components of the radiation reaction potential (Eq. 45).
- `∂Virr_∂a::AbstractArray`: spatial radiation reaction potential derivatives with respect to the harmonic spatial coordinates.
- `∂K_∂xk::AbstractVector{Float64}`: partial derivative of "Kerr potential" K with respect to the harmonic spatial coordinates (Eqs. 54-56, A12-A14).
- `∂Ki_∂xk::AbstractArray`: partial derivative of "Kerr potential" K_i with respect to the harmonic spatial coordinates (Eqs. 54-56, A12-A14).
- `∂Kij_∂xk::AbstractArray`: partial derivative of "Kerr potential" K_ij with respect to the harmonic spatial coordinates (Eqs. 54-56, A12-A14).
- `Q::Float64`: Kerr potential tt component (Eq. 54).
- `Qi::AbstractVector{Float64}`: Kerr potential ti components (Eq. 55).
- `Qij::AbstractArray`: Kerr potential ij (spatial) components (Eq. 56).
- `aSF_H::AbstractArray`: self-acceleration (Eq. 57) in Harmonic coordinates.
- `aSF_BL::AbstractArray`: self-acceleration in Boyer-Lindquist coordinates.
- `a::Float64`: Kerr black hole spin parameter.
"""

# define some useful functions
otimes(a::AbstractVector{Float64}, b::AbstractVector{Float64}) = [a[i] * b[j] for i=1:size(a, 1), j=1:size(b, 1)]    # tensor product of two vectors
otimes(a::AbstractVector{Float64}) = [a[i] * a[j] for i=1:size(a, 1), j=1:size(a, 1)]    # tensor product of a vector with itself
dot3d(u::AbstractVector{Float64}, v::AbstractVector{Float64}) = u[1] * v[1] + u[2] * v[2] + u[3] * v[3]
norm2_3d(u::AbstractVector{Float64}) = u[1] * u[1] + u[2] * u[2] + u[3] * u[3]
norm_3d(u::AbstractVector{Float64}) = sqrt(norm2_3d(u))
dot4d(u::AbstractVector{Float64}, v::AbstractVector{Float64}) = u[1] * v[1] + u[2] * v[2] + u[3] * v[3] + u[4] * v[4]
norm2_4d(u::AbstractVector{Float64}) = u[1] * u[1] + u[2] * u[2] + u[3] * u[3] + u[4] * u[4]
norm_4d(u::AbstractVector{Float64}) = sqrt(norm2_4d(u))

const ημν = [-1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]    # minkowski metric
const ηij = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]    # spatial part of minkowski metric
δ(x::Int, y::Int)::Int = x == y ? 1 : 0

# define vector and scalar potentials for self-force calculation - underscore denotes covariant indices
K(xH::AbstractArray, a::Float64) = HarmonicCoords.g_tt_H(xH, a) + 1.0                         # outputs K00 (Eq. 54)
K_i(xH::AbstractArray, a::Float64) = HarmonicCoords.g_tr_H(xH, a)                             # outputs Ki vector, i.e., Ki for i ∈ {1, 2, 3} (Eq. 55)
K_ij(xH::AbstractArray, a::Float64) = HarmonicCoords.g_rr_H(xH, a) - ηij                      # outputs Kij matrix (Eq. 56)
K_μν(xH::AbstractArray, a::Float64) = HarmonicCoords.g_μν_H(xH, a) - ημν                      # outputs Kμν matrix
Q(xH::AbstractArray, a::Float64) = HarmonicCoords.gTT_H(xH, a) + 1.0                          # outputs Q^00 (Eq. 54)
Qi(xH::AbstractArray, a::Float64) = HarmonicCoords.gTR_H(xH, a)                               # outputs Q^i vector, i.e., Q^i for i ∈ {1, 2, 3} (Eq. 55)
Qij(xH::AbstractArray, a::Float64) = HarmonicCoords.gRR_H(xH, a) - ηij                        # outputs diagonal of Q^ij matrix (Eq. 56)
Qμν(xH::AbstractArray, a::Float64) = HarmonicCoords.gμν_H(xH, a) - ημν                        # outputs Qμν matrix

# define partial derivatives of K (in harmonic coordinates)
# ∂ₖK: outputs float
function ∂K_∂xk(xH::AbstractArray, xBL::AbstractArray, jBLH::AbstractArray, HessBLH::AbstractArray, a::Float64, k::Int)   # Eq. A12
    ∂K=0.0
    @inbounds for μ=1:4
        for i=1:3
            ∂K += 2 * g_μν(xBL[1], xBL[2], xBL[3], a, 1, μ) * Γαμν(xBL[1], xBL[2], xBL[3], a, μ, 1, i+1) * jBLH[i, k]   # i → i + 1 to go from spatial indices to spacetime indices
        end
    end
    return ∂K
end

# ∂ₖKᵢ: outputs float.
function ∂Ki_∂xk(xH::AbstractArray, rH::Float64, xBL::AbstractArray, jBLH::AbstractArray, HessBLH::AbstractArray, a::Float64, k::Int, i::Int)   # Eq. A13
    ∂K=0.0
    @inbounds for m=1:3   # start with iteration over m to not over-count last terms
        ∂K += g_μν(xBL[1], xBL[2], xBL[3], a, 1, m+1) * HessBLH[m][k, i]   # last term Eq. A13, m → m + 1 to go from spatial indices to spacetime indices
        @inbounds for μ=1:4, n=1:3
            ∂K += (g_μν(xBL[1], xBL[2], xBL[3], a, μ, 1) * Γαμν(xBL[1], xBL[2], xBL[3], a, μ, m+1, n+1) + g_μν(xBL[1], xBL[2], xBL[3], a, μ, m+1) * Γαμν(xBL[1], xBL[2], xBL[3], a, μ, 1, n+1)) * jBLH[n, k] * jBLH[m, i]   # first term of Eq. A13
        end
    end
    return ∂K
end

# ∂ₖKᵢⱼ: outputs float.
function ∂Kij_∂xk(xH::AbstractArray, rH::Float64, xBL::AbstractArray, jBLH::AbstractArray, HessBLH::AbstractArray, a::Float64, k::Int, i::Int, j::Int)   # Eq. A14
    ∂K=0.0
    @inbounds for m=1:3
        @inbounds for l=1:3   # iterate over m and l first to avoid over-counting
            ∂K += g_μν(xBL[1], xBL[2], xBL[3], a, l+1, m+1) * (HessBLH[l][k, i] * jBLH[m, j] + HessBLH[l][k, j] * jBLH[m, i])  # last term Eq. A14
            @inbounds for μ=1:4, n=1:3
                ∂K += (g_μν(xBL[1], xBL[2], xBL[3], a, μ, l+1) * Γαμν(xBL[1], xBL[2], xBL[3], a, μ, m+1, n+1) + g_μν(xBL[1], xBL[2], xBL[3], a, μ, m+1) * Γαμν(xBL[1], xBL[2], xBL[3], a, μ, l+1, n+1)) * jBLH[n, k] * jBLH[l, i] * jBLH[m, j]   # first term of Eq. A14
            end
        end
    end
    return ∂K
end

# define relativistic Γ factor
Γ(vH::AbstractArray, xH::AbstractArray, a::Float64) = 1.0 / sqrt(1.0 - SelfAccelerationHarmonic2.norm2_3d(vH) - K(xH, a) - 2.0 * dot(K_i(xH, a), vH) - transpose(vH) * K_ij(xH, a) * vH)   # Eq. A3

# define projection operator
Pαβ(vH::AbstractArray, xH::AbstractArray, a::Float64) = ημν + Qμν(xH, a) + Γ(vH, xH, a)^2 * otimes(vcat([1], vH))   # contravariant, Eq. A1
P_αβ(vH::AbstractArray, xH::AbstractArray, a::Float64) =  ημν + K_μν(xH, a) + Γ(vH, xH, a)^2 * otimes(vcat([-1], vH))   # cοvariant, Eq. A2 (note that we take both contravariant and covariant velocities as arguments)

### SELF-ACCELERATION PIECES ###
# compute self-acceleration pieces. crucially ∂Virr_∂a = ∂_{a}V_{i}
function A_RR(v::Float64, vH::AbstractArray, ∂Vrr_∂t::Float64, ∂Vrr_∂a::MVector{3, Float64}, ∂Virr_∂a::MMatrix{3, 3, Float64})
    aRR = (1.0 - v^2) * ∂Vrr_∂t   # first term in Eq. A4
    @inbounds for i=1:3
        aRR += 2.0 * vH[i] * ∂Vrr_∂a[i]   # second term Eq. A4
        @inbounds for j=1:3
            aRR += -4.0 * vH[i] * vH[j] * ∂Virr_∂a[j, i]   # third term Eq. A4
        end
    end
    return aRR
end

function Ai_RR(v::Float64,  vH::AbstractArray, ∂Vrr_∂t::Float64, ∂Virr_∂t::MVector{3, Float64}, ∂Vrr_∂a::MVector{3, Float64}, ∂Virr_∂a::MMatrix{3, 3, Float64}, i::Int)
    aiRR = -(1 + v^2) * ∂Vrr_∂a[i] + 2.0 * vH[i] * ∂Vrr_∂t - 4.0 * ∂Virr_∂t[i]    # first, second, and last term in Eq. A5
    @inbounds for j=1:3
        aiRR += 2.0 * vH[i] * vH[j] * ∂Vrr_∂a[j] - 4.0 * vH[j] * (∂Virr_∂a[i, j] - ∂Virr_∂a[j, i])    # third and fourth terms in Eq. A5
    end
    return aiRR
end

function A1_β(v::Float64, vH::AbstractArray, ∂Vrr_∂t::Float64, ∂Vrr_∂a::MVector{3, Float64}, ∂Virr_∂t::MVector{3, Float64}, ∂Virr_∂a::MMatrix{3, 3, Float64, 9})
    return [i==1 ? A_RR(v, vH, ∂Vrr_∂t, ∂Vrr_∂a, ∂Virr_∂a) : Ai_RR(v, vH, ∂Vrr_∂t, ∂Virr_∂t, ∂Vrr_∂a, ∂Virr_∂a, i-1) for i = 1:4]
end

function B_RR(Qi::AbstractVector{Float64}, ∂K_∂xk::SVector{3, Float64})
    return dot(Qi, ∂K_∂xk)   # Eq. A6
end

function Bi_RR(Qij::AbstractArray, ∂K_∂xk::SVector{3, Float64})
    return -2.0 * (ηij + Qij) * ∂K_∂xk   # Eq. A9
end

# Eq. A7
function C_RR(vH::AbstractArray, ∂K_∂xk::SVector{3, Float64}, ∂Ki_∂xk::SMatrix{3, 3, Float64}, Q::Float64, Qi::AbstractVector{Float64})
    C = 0.0
    @inbounds for i=1:3
        C += 2.0 * (1.0 - Q) * vH[i] * ∂K_∂xk[i]
        @inbounds for j=1:3
            C += 2.0 * Qi[i] * vH[j] * (∂Ki_∂xk[i, j] - ∂Ki_∂xk[j, i])
        end
    end
    return C
end

# Eq. A10
function Ci_RR(vH::AbstractArray, ∂K_∂xk::SVector{3, Float64}, ∂Ki_∂xk::SMatrix{3, 3, Float64}, Qi::AbstractVector{Float64}, Qij::AbstractArray)   # Eq. A10
    C = @MVector [0., 0., 0.]
    @inbounds for j=1:3
        @inbounds for i=1:3
            C[i] += 4.0 * Qi[i] * vH[j] * ∂K_∂xk[j]
        end
        C .+= 4.0 * (ηij + Qij) * vH[j] * ([(∂Ki_∂xk[j, k] - ∂Ki_∂xk[k, j]) for k=1:3]) 
    end
    return C
end

# Eq. A8
function D_RR(vH::AbstractArray, ∂Ki_∂xk::SMatrix{3, 3, Float64}, ∂Kij_∂xk::SArray{Tuple{3, 3, 3}, Float64, 3, 27}, Q::Float64, Qi::AbstractVector{Float64})
    D = 0.0
    @inbounds for i=1:3
        @inbounds for j=1:3
            D += 2.0 * (1.0 - Q) * vH[i] * vH[j] * ∂Ki_∂xk[i, j]
            @inbounds for k=1:3
                D += -Qi[i] * vH[j] * vH[k] * (∂Kij_∂xk[j, k, i] + ∂Kij_∂xk[k, j, i] - ∂Kij_∂xk[i, j, k]) 
            end
        end
    end
    return D
end

# Eq. A11
function Di_RR(vH::AbstractArray, ∂Ki_∂xk::SMatrix{3, 3, Float64}, ∂Kij_∂xk::SArray{Tuple{3, 3, 3}, Float64, 3, 27}, Qi::AbstractVector{Float64}, Qij::AbstractArray)   # Eq. A11
    D = @MVector [0., 0., 0.]
    @inbounds for j=1:3
        @inbounds for k=1:3
            @inbounds for i=1:3
                D[i] += 4.0 * Qi[i] * vH[j] * vH[k] * ∂Ki_∂xk[j, k]
            end
            D .+= 2.0 * (ηij + Qij) * vH[j] * vH[k] * [(∂Kij_∂xk[j, k, l] + ∂Kij_∂xk[k, j, l] - ∂Kij_∂xk[l, j, k]) for l=1:3]
        end
    end
    return D
end


# computes the four self-acceleration components A^{2}_{β} (Eqs. 62 - 63)
function A2_β(xH::AbstractArray, vH::AbstractArray, xBL::AbstractArray, rH::Float64, a::Float64, VRR::Float64, ViRR::MVector{3, Float64})
    jBLH = HarmonicCoords.jBLH(xH, a)
    HessBLH = [HarmonicCoords.HessBLH(xH, rH, a, m) for m=1:3]
    ∂K_∂xk = @SVector [SelfAccelerationHarmonic2.∂K_∂xk(xH, xBL, jBLH, HessBLH, a, j) for j=1:3];
    ∂Ki_∂xk = @SMatrix [SelfAccelerationHarmonic2.∂Ki_∂xk(xH, rH, xBL, jBLH, HessBLH, a, j, k) for j=1:3, k=1:3];
    ∂Kij_∂xk = @SArray [SelfAccelerationHarmonic2.∂Kij_∂xk(xH, rH, xBL, jBLH, HessBLH, a, j, k, l) for j=1:3, k=1:3, l=1:3]
    Q = SelfAccelerationHarmonic2.Q(xH, a)
    Qi = SelfAccelerationHarmonic2.Qi(xH, a)
    Qij = SelfAccelerationHarmonic2.Qij(xH, a)

    BRR = B_RR(Qi, ∂K_∂xk)
    BiRR = Bi_RR(Qij, ∂K_∂xk)

    CRR = C_RR(vH, ∂K_∂xk, ∂Ki_∂xk, Q, Qi)
    CiRR = Ci_RR(vH, ∂K_∂xk, ∂Ki_∂xk, Qi, Qij)

    DRR = D_RR(vH, ∂Ki_∂xk, ∂Kij_∂xk, Q, Qi)
    DiRR = Di_RR(vH, ∂Ki_∂xk, ∂Kij_∂xk, Qi, Qij)

    A2_t = (BRR + CRR + DRR) * VRR + dot((BiRR + CiRR + DiRR), ViRR)   # Eq. 62
    A2_i = -2.0 * (BRR + CRR + DRR) * ViRR - (BiRR + CiRR + DiRR) * VRR / 2.0  # Eq. 63

    return vcat(A2_t, A2_i)
end

# compute self-acceleration in harmonic coordinates and transform components back to BL
function aRRα(aSF_H::AbstractVector{Float64}, aSF_BL::AbstractVector{Float64}, xH::AbstractVector{Float64}, v::Float64, vH::AbstractVector{Float64}, xBL::AbstractVector{Float64}, rH::Float64, a::Float64, Vrr::Float64, ∂Vrr_∂t::Float64, Virr::MVector{3, Float64}, ∂Vrr_∂a::MVector{3, Float64}, ∂Virr_∂t::MVector{3, Float64}, ∂Virr_∂a::MMatrix{3, 3, Float64, 9})
    aSF_H[:] = -Γ(vH, xH, a)^2 * HarmonicCoords.gμν_H(xH, a) * (A1_β(v, vH, ∂Vrr_∂t, ∂Vrr_∂a, ∂Virr_∂t, ∂Virr_∂a) - A2_β(xH, vH, xBL, rH, a, Vrr, Virr))
    aSF_BL[1] = aSF_H[1]
    aSF_BL[2:4] = HarmonicCoords.aHtoBL(xH, zeros(3), aSF_H[2:4], a)
end

end