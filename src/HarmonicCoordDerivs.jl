module HarmonicCoordDerivs

const M::Float64 = 1.0

ρ(r::Float64, a::Float64) = sqrt(a^2 + (-M + r)^2)

d1ρ(r::Float64, a::Float64) = (-M + r)/sqrt(a^2 + (M - r)^2)

d2ρ(r::Float64, a::Float64) = a^2/(a^2 + (M - r)^2)^1.5

d3ρ(r::Float64, a::Float64) = (3*a^2*(M - r))/(a^2 + (M - r)^2)^2.5

d4ρ(r::Float64, a::Float64) = (-3*a^2*(a^2 - 4*(M - r)^2))/(a^2 + (M - r)^2)^3.5

d5ρ(r::Float64, a::Float64) = (15*a^2*(-3*a^2 + 4*(M - r)^2)*(M - r))/(a^2 + (M - r)^2)^4.5

d6ρ(r::Float64, a::Float64) = (45*a^2*(a^4 - 12*a^2*(M - r)^2 + 8*(M - r)^4))/(a^2 + (M - r)^2)^5.5

d7ρ(r::Float64, a::Float64) = (315*a^2*(5*a^4 - 20*a^2*(M - r)^2 + 8*(M - r)^4)*(M - r))/(a^2 + (M - r)^2)^6.5

d8ρ(r::Float64, a::Float64) = (315*a^2*(-5*a^6 + 120*a^4*(M - r)^2 - 240*a^2*(M - r)^4 + 64*(M - r)^6))/(a^2 + (M - r)^2)^7.5

d9ρ(r::Float64, a::Float64) = (2835*a^2*(-35*a^6 + 280*a^4*(M - r)^2 - 336*a^2*(M - r)^4 + 64*(M - r)^6)*(M - r))/(a^2 + (M - r)^2)^8.5

function compute_ρ_derivs!(dρ::AbstractVector{Float64}, r::Float64, a::Float64)
    dρ[1] = d1ρ(r, a)
    dρ[2] = d2ρ(r, a)
    dρ[3] = d3ρ(r, a)
    dρ[4] = d4ρ(r, a)
    dρ[5] = d5ρ(r, a)
    dρ[6] = d6ρ(r, a)
    dρ[7] = d7ρ(r, a)
    dρ[8] = d8ρ(r, a)
    dρ[9] = d9ρ(r, a)
end

Φ(r::Float64, a::Float64, rm::Float64, rp::Float64) = π/2. + atan(M - r, a) + (a*log((r - rp)/(r - rm)))/(2.0*sqrt(-a^2 + M^2))

d1Φ(r::Float64, a::Float64, rm::Float64, rp::Float64) = -(a/(a^2 + (M - r)^2)) + (a*(-rm + rp))/(2.0*sqrt(-a^2 + M^2)*(r - rm)*(r - rp))

d2Φ(r::Float64, a::Float64, rm::Float64, rp::Float64) = (2*a*(-M + r))/(a^2 + (M - r)^2)^2 + (a*(2*r - rm - rp)*(rm - rp))/(2.0*sqrt(-a^2 + M^2)*(r - rm)^2*(r - rp)^2)

d3Φ(r::Float64, a::Float64, rm::Float64, rp::Float64) = (2*a*(a^2 - 3*(M - r)^2))/(a^2 + (M - r)^2)^3 - (a*(rm - rp)*(3*r^2 + rm^2 + rm*rp + rp^2 - 3*r*(rm + rp)))/(sqrt(-a^2 + M^2)*(r - rm)^3*(r - rp)^3)

d4Φ(r::Float64, a::Float64, rm::Float64, rp::Float64) = (24*a*(M - r)*(a + M - r)*(a - M + r))/(a^2 + (M - r)^2)^4 + (3*a*(2*r - rm - rp)*(rm - rp)*(2*r^2 + rm^2 + rp^2 - 2*r*(rm + rp)))/(sqrt(-a^2 + M^2)*(r - rm)^4*(r - rp)^4)

d5Φ(r::Float64, a::Float64, rm::Float64, rp::Float64) = (-24*a*(a^4 - 10*a^2*(M - r)^2 + 5*(M - r)^4))/(a^2 + (M - r)^2)^5 - (12*a*(rm^5 + 5*r^4*(rm - rp) - rp^5 + 10*r^3*(-rm^2 + rp^2) + 10*r^2*(rm^3 - rp^3) + 5*r*(-rm^4 + rp^4)))/(sqrt(-a^2 + M^2)*(r - rm)^5*(r - rp)^5)

d6Φ(r::Float64, a::Float64, rm::Float64, rp::Float64) = (-240*a*(3*a^4 - 10*a^2*(M - r)^2 + 3*(M - r)^4)*(M - r))/(a^2 + (M - r)^2)^6 + (60*a*(2*r - rm - rp)*(rm - rp)*(3*r^2 + rm^2 + rm*rp + rp^2 - 3*r*(rm + rp))*(r^2 + rm^2 - rm*rp + rp^2 - r*(rm + rp)))/(sqrt(-a^2 + M^2)*(r - rm)^6*(r - rp)^6)

d7Φ(r::Float64, a::Float64, rm::Float64, rp::Float64) = (720*a*(a^6 - 21*a^4*(M - r)^2 + 35*a^2*(M - r)^4 - 7*(M - r)^6))/(a^2 + (M - r)^2)^7 - (360*a*(rm^7 + 7*r^6*(rm - rp) - rp^7 + 21*r^5*(-rm^2 + rp^2) + 35*r^4*(rm^3 - rp^3) + 35*r^3*(-rm^4 + rp^4) + 21*r^2*(rm^5 - rp^5) + 7*r*(-rm^6 + rp^6)))/(sqrt(-a^2 + M^2)*(r - rm)^7*(r - rp)^7)

d8Φ(r::Float64, a::Float64, rm::Float64, rp::Float64) = (40320*a*(a^6 - 7*a^4*(M - r)^2 + 7*a^2*(M - r)^4 - (M - r)^6)*(M - r))/(a^2 + (M - r)^2)^8 + (2520*a*(2*r - rm - rp)*(rm - rp)*(2*r^2 + rm^2 + rp^2 - 2*r*(rm + rp))*(2*r^4 + rm^4 + rp^4 - 4*r^3*(rm + rp) + 6*r^2*(rm^2 + rp^2) - 4*r*(rm^3 + rp^3)))/(sqrt(-a^2 + M^2)*(r - rm)^8*(r - rp)^8)

d9Φ(r::Float64, a::Float64, rm::Float64, rp::Float64) = (-40320*a*(a^8 - 36*a^6*(M - r)^2 + 126*a^4*(M - r)^4 - 84*a^2*(M - r)^6 + 9*(M - r)^8))/(a^2 + (M - r)^2)^9 - (20160*a*(rm^9 + 9*r^8*(rm - rp) - rp^9 + 36*r^7*(-rm^2 + rp^2) + 84*r^6*(rm^3 - rp^3) + 126*r^5*(-rm^4 + rp^4) + 126*r^4*(rm^5 - rp^5) + 84*r^3*(-rm^6 + rp^6) + 36*r^2*(rm^7 - rp^7) + 9*r*(-rm^8 + rp^8)))/(sqrt(-a^2 + M^2)*(r - rm)^9*(r - rp)^9)

# Φ(r::Float64, a::Float64, rm::Float64, rp::Float64) = π/2 - atan(-M + r, a) - (a*log((r - rm)/(r - rp)))/(2.0*sqrt(-a^2 + M^2))

# d1Φ(r::Float64, a::Float64, rm::Float64, rp::Float64) = -(a/(a^2 + (M - r)^2)) + (a*(-rm + rp))/(2.0*sqrt(-a^2 + M^2)*(r - rm)*(r - rp))

# d2Φ(r::Float64, a::Float64, rm::Float64, rp::Float64) = (2*a*(-M + r))/(a^2 + (M - r)^2)^2 + (a*(2*r - rm - rp)*(rm - rp))/(2.0*sqrt(-a^2 + M^2)*(r -
# rm)^2*(r - rp)^2)

# d3Φ(r::Float64, a::Float64, rm::Float64, rp::Float64) = (2*a*(a^2 - 3*(M - r)^2))/(a^2 + (M - r)^2)^3 - (a*(rm - rp)*(3*r^2 + rm^2 + rm*rp + rp^2 -
# 3*r*(rm + rp)))/(sqrt(-a^2 + M^2)*(r - rm)^3*(r - rp)^3)

# d4Φ(r::Float64, a::Float64, rm::Float64, rp::Float64) = (24*a*(M - r)*(a + M - r)*(a - M + r))/(a^2 + (M - r)^2)^4 + (3*a*(2*r - rm - rp)*(rm - rp)*(2*r^2 +
# rm^2 + rp^2 - 2*r*(rm + rp)))/(sqrt(-a^2 + M^2)*(r - rm)^4*(r - rp)^4)

# d5Φ(r::Float64, a::Float64, rm::Float64, rp::Float64) = (-24*a*(a^4 - 10*a^2*(M - r)^2 + 5*(M - r)^4))/(a^2 + (M - r)^2)^5 - (12*a*(rm^5 + 5*r^4*(rm - rp) -
# rp^5 + 10*r^3*(-rm^2 + rp^2) + 10*r^2*(rm^3 - rp^3) + 5*r*(-rm^4 + rp^4)))/(sqrt(-a^2 + M^2)*(r - rm)^5*(r - rp)^5)

# d6Φ(r::Float64, a::Float64, rm::Float64, rp::Float64) = (-240*a*(3*a^4 - 10*a^2*(M - r)^2 + 3*(M - r)^4)*(M - r))/(a^2 + (M - r)^2)^6 +
# (60*a*(2*r - rm - rp)*(rm - rp)*(3*r^2 + rm^2 + rm*rp + rp^2 - 3*r*(rm + rp))*(r^2 + rm^2 - rm*rp + rp^2 - r*(rm + rp)))/(sqrt(-a^2 + M^2)*(r - rm)^6*(r - rp)^6)

# d7Φ(r::Float64, a::Float64, rm::Float64, rp::Float64) = (720*a*(a^6 - 21*a^4*(M - r)^2 + 35*a^2*(M - r)^4 - 7*(M - r)^6))/(a^2 + (M - r)^2)^7 -
# (360*a*(rm^7 + 7*r^6*(rm - rp) - rp^7 + 21*r^5*(-rm^2 + rp^2) + 35*r^4*(rm^3 - rp^3) + 35*r^3*(-rm^4 + rp^4) + 21*r^2*(rm^5 - rp^5) + 7*r*(-rm^6 + rp^6)))/(sqrt(-a^2 + M^2)*(r - rm)^7*(r - rp)^7)

# d8Φ(r::Float64, a::Float64, rm::Float64, rp::Float64) = (40320*a*(a^6 - 7*a^4*(M - r)^2 + 7*a^2*(M - r)^4 - (M - r)^6)*(M - r))/(a^2 + (M - r)^2)^8 +
# (2520*a*(2*r - rm - rp)*(rm - rp)*(2*r^2 + rm^2 + rp^2 - 2*r*(rm + rp))*(2*r^4 + rm^4 + rp^4 - 4*r^3*(rm + rp) + 6*r^2*(rm^2 + rp^2) -
# 4*r*(rm^3 + rp^3)))/(sqrt(-a^2 + M^2)*(r - rm)^8*(r - rp)^8)

# d9Φ(r::Float64, a::Float64, rm::Float64, rp::Float64) = (-40320*a*(a^8 - 36*a^6*(M - r)^2 + 126*a^4*(M - r)^4 - 84*a^2*(M - r)^6 + 9*(M - r)^8))/(a^2 + (M - r)^2)^9 - (20160*a*(rm^9 + 9*r^8*(rm - rp) - rp^9 + 36*r^7*(-rm^2 + rp^2) + 84*r^6*(rm^3 - rp^3) + 126*r^5*(-rm^4 + rp^4) + 126*r^4*(rm^5 - rp^5) + 84*r^3*(-rm^6 + rp^6) + 36*r^2*(rm^7 - rp^7) + 9*r*(-rm^8 + rp^8)))/(sqrt(-a^2 + M^2)*(r - rm)^9*(r - rp)^9)

function compute_Φ_derivs!(dΦ::AbstractVector{Float64}, r::Float64, a::Float64, rm::Float64, rp::Float64)
    dΦ[1] = d1Φ(r, a, rm, rp)
    dΦ[2] = d2Φ(r, a, rm, rp)
    dΦ[3] = d3Φ(r, a, rm, rp)
    dΦ[4] = d4Φ(r, a, rm, rp)
    dΦ[5] = d5Φ(r, a, rm, rp)
    dΦ[6] = d6Φ(r, a, rm, rp)
    dΦ[7] = d7Φ(r, a, rm, rp)
    dΦ[8] = d8Φ(r, a, rm, rp)
    dΦ[9] = d9Φ(r, a, rm, rp)
end

ξ(x::AbstractVector{Float64}, Φ::Float64) = x[3] - Φ

d1ξ(dx::AbstractVector{Float64}, dΦ::AbstractVector{Float64}) = dx[3] - dx[1]*dΦ[1]

d2ξ(dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, dΦ::AbstractVector{Float64}) = -(dΦ[1]*d2x[1]) + d2x[3] - dx[1]^2*dΦ[2]

d3ξ(dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, dΦ::AbstractVector{Float64}) = -3*dx[1]*d2x[1]*dΦ[2] - dΦ[1]*d3x[1] + d3x[3] - dx[1]^3*dΦ[3]

d4ξ(dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, d4x::AbstractVector{Float64}, dΦ::AbstractVector{Float64}) = -3*d2x[1]^2*dΦ[2] - 4*dx[1]*dΦ[2]*d3x[1] -
6*dx[1]^2*d2x[1]*dΦ[3] - dΦ[1]*d4x[1] + d4x[3] - dx[1]^4*dΦ[4]

d5ξ(dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, d4x::AbstractVector{Float64}, d5x::AbstractVector{Float64}, dΦ::AbstractVector{Float64}) = -10*d2x[1]*dΦ[2]*d3x[1] -
15*dx[1]*d2x[1]^2*dΦ[3] - 10*dx[1]^2*d3x[1]*dΦ[3] - 5*dx[1]*dΦ[2]*d4x[1] - 10*dx[1]^3*d2x[1]*dΦ[4] - dΦ[1]*d5x[1] + d5x[3] - dx[1]^5*dΦ[5]

d6ξ(dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, d4x::AbstractVector{Float64}, d5x::AbstractVector{Float64}, d6x::AbstractVector{Float64}, dΦ::AbstractVector{Float64}) = -10*dΦ[2]*d3x[1]^2 -
15*d2x[1]^3*dΦ[3] - 60*dx[1]*d2x[1]*d3x[1]*dΦ[3] - 15*d2x[1]*dΦ[2]*d4x[1] - 15*dx[1]^2*dΦ[3]*d4x[1] - 45*dx[1]^2*d2x[1]^2*dΦ[4] - 20*dx[1]^3*d3x[1]*dΦ[4] - 6*dx[1]*dΦ[2]*d5x[1] - 15*dx[1]^4*d2x[1]*dΦ[5] - dΦ[1]*d6x[1] + d6x[3] - dx[1]^6*dΦ[6]

d7ξ(dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, d4x::AbstractVector{Float64}, d5x::AbstractVector{Float64}, d6x::AbstractVector{Float64}, d7x::AbstractVector{Float64},
dΦ::AbstractVector{Float64}) = -105*d2x[1]^2*d3x[1]*dΦ[3] - 70*dx[1]*d3x[1]^2*dΦ[3] - 35*dΦ[2]*d3x[1]*d4x[1] - 105*dx[1]*d2x[1]*dΦ[3]*d4x[1] - 105*dx[1]*d2x[1]^3*dΦ[4] -
210*dx[1]^2*d2x[1]*d3x[1]*dΦ[4] - 35*dx[1]^3*d4x[1]*dΦ[4] - 21*d2x[1]*dΦ[2]*d5x[1] - 21*dx[1]^2*dΦ[3]*d5x[1] - 105*dx[1]^3*d2x[1]^2*dΦ[5] - 35*dx[1]^4*d3x[1]*dΦ[5] -
7*dx[1]*dΦ[2]*d6x[1] - 21*dx[1]^5*d2x[1]*dΦ[6] - dΦ[1]*d7x[1] + d7x[3] - dx[1]^7*dΦ[7]

d8ξ(dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, d4x::AbstractVector{Float64}, d5x::AbstractVector{Float64}, d6x::AbstractVector{Float64}, d7x::AbstractVector{Float64}, d8x::AbstractVector{Float64},
dΦ::AbstractVector{Float64}) = -280*d2x[1]*d3x[1]^2*dΦ[3] - 210*d2x[1]^2*dΦ[3]*d4x[1] - 280*dx[1]*d3x[1]*dΦ[3]*d4x[1] - 35*dΦ[2]*d4x[1]^2 - 105*d2x[1]^4*dΦ[4] -
840*dx[1]*d2x[1]^2*d3x[1]*dΦ[4] - 280*dx[1]^2*d3x[1]^2*dΦ[4] - 420*dx[1]^2*d2x[1]*d4x[1]*dΦ[4] - 56*dΦ[2]*d3x[1]*d5x[1] - 168*dx[1]*d2x[1]*dΦ[3]*d5x[1] -
56*dx[1]^3*dΦ[4]*d5x[1] - 420*dx[1]^2*d2x[1]^3*dΦ[5] - 560*dx[1]^3*d2x[1]*d3x[1]*dΦ[5] - 70*dx[1]^4*d4x[1]*dΦ[5] - 28*d2x[1]*dΦ[2]*d6x[1] - 28*dx[1]^2*dΦ[3]*d6x[1] -
210*dx[1]^4*d2x[1]^2*dΦ[6] - 56*dx[1]^5*d3x[1]*dΦ[6] - 8*dx[1]*dΦ[2]*d7x[1] - 28*dx[1]^6*d2x[1]*dΦ[7] - dΦ[1]*d8x[1] + d8x[3] - dx[1]^8*dΦ[8]

# d9ξ(dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, d4x::AbstractVector{Float64}, d5x::AbstractVector{Float64}, d6x::AbstractVector{Float64}, d7x::AbstractVector{Float64}, d8x::AbstractVector{Float64}, d9x::AbstractVector{Float64}, dΦ::AbstractVector{Float64}) = -280*d3x[1]^3*dΦ[3] - 1260*d2x[1]*d3x[1]*dΦ[3]*d4x[1] - 315*dx[1]*dΦ[3]*d4x[1]^2 - 1260*d2x[1]^3*d3x[1]*dΦ[4] - 2520*dx[1]*d2x[1]*d3x[1]^2*dΦ[4] - 1890*dx[1]*d2x[1]^2*d4x[1]*dΦ[4] - 1260*dx[1]^2*d3x[1]*d4x[1]*dΦ[4] - 378*d2x[1]^2*dΦ[3]*d5x[1] - 504*dx[1]*d3x[1]*dΦ[3]*d5x[1] - 126*dΦ[2]*d4x[1]*d5x[1] - 756*dx[1]^2*d2x[1]*dΦ[4]*d5x[1] - 945*dx[1]*d2x[1]^4*dΦ[5] - 3780*dx[1]^2*d2x[1]^2*d3x[1]*dΦ[5] - 840*dx[1]^3*d3x[1]^2*dΦ[5] - 1260*dx[1]^3*d2x[1]*d4x[1]*dΦ[5] - 126*dx[1]^4*d5x[1]*dΦ[5] - 84*dΦ[2]*d3x[1]*d6x[1] - 252*dx[1]*d2x[1]*dΦ[3]*d6x[1] - 84*dx[1]^3*dΦ[4]*d6x[1] - 1260*dx[1]^3*d2x[1]^3*dΦ[6] - 1260*dx[1]^4*d2x[1]*d3x[1]*dΦ[6] - 126*dx[1]^5*d4x[1]*dΦ[6] - 36*d2x[1]*dΦ[2]*d7x[1] - 36*dx[1]^2*dΦ[3]*d7x[1] - 378*dx[1]^5*d2x[1]^2*dΦ[7] - 84*dx[1]^6*d3x[1]*dΦ[7] - 9*dx[1]*dΦ[2]*d8x[1] - 36*dx[1]^7*d2x[1]*dΦ[8] - dΦ[1]*d9x[1] + dΦ[9] - dx[1]^9*dΦ[9]

d9ξ(dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, d4x::AbstractVector{Float64}, d5x::AbstractVector{Float64}, d6x::AbstractVector{Float64}, d7x::AbstractVector{Float64}, d8x::AbstractVector{Float64}, d9x::AbstractVector{Float64}, dΦ::AbstractVector{Float64}) = -280*d3x[1]^3*dΦ[3] - 1260*d2x[1]*d3x[1]*dΦ[3]*d4x[1] - 315*dx[1]*dΦ[3]*d4x[1]^2 - 1260*d2x[1]^3*d3x[1]*dΦ[4] - 2520*dx[1]*d2x[1]*d3x[1]^2*dΦ[4] - 1890*dx[1]*d2x[1]^2*d4x[1]*dΦ[4] - 1260*dx[1]^2*d3x[1]*d4x[1]*dΦ[4] - 378*d2x[1]^2*dΦ[3]*d5x[1] - 504*dx[1]*d3x[1]*dΦ[3]*d5x[1] - 126*dΦ[2]*d4x[1]*d5x[1] - 756*dx[1]^2*d2x[1]*dΦ[4]*d5x[1] - 945*dx[1]*d2x[1]^4*dΦ[5] - 3780*dx[1]^2*d2x[1]^2*d3x[1]*dΦ[5] - 840*dx[1]^3*d3x[1]^2*dΦ[5] - 1260*dx[1]^3*d2x[1]*d4x[1]*dΦ[5] - 126*dx[1]^4*d5x[1]*dΦ[5] - 84*dΦ[2]*d3x[1]*d6x[1] - 252*dx[1]*d2x[1]*dΦ[3]*d6x[1] - 84*dx[1]^3*dΦ[4]*d6x[1] - 1260*dx[1]^3*d2x[1]^3*dΦ[6] - 1260*dx[1]^4*d2x[1]*d3x[1]*dΦ[6] - 126*dx[1]^5*d4x[1]*dΦ[6] - 36*d2x[1]*dΦ[2]*d7x[1] - 36*dx[1]^2*dΦ[3]*d7x[1] - 378*dx[1]^5*d2x[1]^2*dΦ[7] - 84*dx[1]^6*d3x[1]*dΦ[7] - 9*dx[1]*dΦ[2]*d8x[1] - 36*dx[1]^7*d2x[1]*dΦ[8] - dΦ[1]*d9x[1] + d9x[3] - dx[1]^9*dΦ[9]

function compute_ξ_derivs!(dξ::AbstractVector{Float64}, dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, d4x::AbstractVector{Float64}, d5x::AbstractVector{Float64}, d6x::AbstractVector{Float64}, d7x::AbstractVector{Float64}, d8x::AbstractVector{Float64}, d9x::AbstractVector{Float64}, dΦ::AbstractVector{Float64})
    dξ[1] = d1ξ(dx, dΦ)
    dξ[2] = d2ξ(dx, d2x, dΦ)
    dξ[3] = d3ξ(dx, d2x, d3x, dΦ)
    dξ[4] = d4ξ(dx, d2x, d3x, d4x, dΦ)
    dξ[5] = d5ξ(dx, d2x, d3x, d4x, d5x, dΦ)
    dξ[6] = d6ξ(dx, d2x, d3x, d4x, d5x, d6x, dΦ)
    dξ[7] = d7ξ(dx, d2x, d3x, d4x, d5x, d6x, d7x, dΦ)
    dξ[8] = d8ξ(dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, dΦ)
    dξ[9] = d9ξ(dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, dΦ)
end

xH(x::AbstractVector{Float64}, ξ::Float64, ρ::Float64) = sin(x[2])*cos(ξ)*ρ

dxH1(x::AbstractVector{Float64}, dx::AbstractVector{Float64}, ξ::Float64, dξ::AbstractVector{Float64}, ρ::Float64, dρ::AbstractVector{Float64}) = cos(x[2])*cos(ξ)*ρ*dx[2] - sin(x[2])*sin(ξ)*ρ*dξ[1] +
cos(ξ)*sin(x[2])*dx[1]*dρ[1]

dxH2(x::AbstractVector{Float64}, dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, ξ::Float64, dξ::AbstractVector{Float64}, ρ::Float64, dρ::AbstractVector{Float64}) = 2*dx[1]*(cos(x[2])*cos(ξ)*dx[2] -
sin(x[2])*sin(ξ)*dξ[1])*dρ[1] + ρ*(-2*cos(x[2])*sin(ξ)*dx[2]*dξ[1] + cos(ξ)*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) + sin(x[2])*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2])) +
cos(ξ)*sin(x[2])*(dρ[1]*d2x[1] + dx[1]^2*dρ[2])

dxH3(x::AbstractVector{Float64}, dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, ξ::Float64, dξ::AbstractVector{Float64}, ρ::Float64, dρ::AbstractVector{Float64}) = 3*dx[1]*dρ[1]*(-2*
cos(x[2])*sin(ξ)*dx[2]*dξ[1] + cos(ξ)*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) + sin(x[2])*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2])) + 3*(cos(x[2])*cos(ξ)*dx[2] -
sin(x[2])*sin(ξ)*dξ[1])*(dρ[1]*d2x[1] + dx[1]^2*dρ[2]) + ρ*(-3*sin(ξ)*dξ[1]*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) + 3*cos(x[2])*dx[2]*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2]) +
cos(ξ)*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) + sin(x[2])*(sin(ξ)*dξ[1]^3 - 3*cos(ξ)*dξ[1]*dξ[2] - sin(ξ)*dξ[3])) + cos(ξ)*sin(x[2])*(3*dx[1]*
d2x[1]*dρ[2] + dρ[1]*d3x[1] + dx[1]^3*dρ[3])

dxH4(x::AbstractVector{Float64}, dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, d4x::AbstractVector{Float64}, ξ::Float64, dξ::AbstractVector{Float64}, ρ::Float64, dρ::AbstractVector{Float64}) = 6*(-2*cos(x[2])*sin(ξ)*dx[2]*dξ[1] + cos(ξ)*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) + sin(x[2])*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2]))*(dρ[1]*d2x[1] + dx[1]^2*dρ[2]) +
4*dx[1]*dρ[1]*(-3*sin(ξ)*dξ[1]*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) + 3*cos(x[2])*dx[2]*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2]) + cos(ξ)*(-(cos(x[2])*dx[2]^3) -
3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) + sin(x[2])*(sin(ξ)*dξ[1]^3 - 3*cos(ξ)*dξ[1]*dξ[2] - sin(ξ)*dξ[3])) + 4*(cos(x[2])*cos(ξ)*dx[2] -
sin(x[2])*sin(ξ)*dξ[1])*(3*dx[1]*d2x[1]*dρ[2] + dρ[1]*d3x[1] + dx[1]^3*dρ[3]) + ρ*(6*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2]) -
4*sin(ξ)*dξ[1]*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) + 4*cos(x[2])*dx[2]*(sin(ξ)*dξ[1]^3 - 3*cos(ξ)*dξ[1]*dξ[2] - sin(ξ)*dξ[3]) +
cos(ξ)*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + sin(x[2])*(cos(ξ)*dξ[1]^4 +
6*sin(ξ)*dξ[1]^2*dξ[2] - 3*cos(ξ)*dξ[2]^2 - 4*cos(ξ)*dξ[1]*dξ[3] - sin(ξ)*dξ[4])) + cos(ξ)*sin(x[2])*(3*d2x[1]^2*dρ[2] + 4*dx[1]*dρ[2]*d3x[1] + 6*dx[1]^2*d2x[1]*dρ[3] +
dρ[1]*d4x[1] + dx[1]^4*dρ[4])

dxH5(x::AbstractVector{Float64}, dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, d4x::AbstractVector{Float64}, d5x::AbstractVector{Float64}, ξ::Float64, dξ::AbstractVector{Float64}, ρ::Float64, dρ::AbstractVector{Float64}) = 10*(dρ[1]*d2x[1] + dx[1]^2*dρ[2])*(-3*sin(ξ)*dξ[1]*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) + 3*cos(x[2])*dx[2]*(-(cos(ξ)*dξ[1]^2) -
sin(ξ)*dξ[2]) + cos(ξ)*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) + sin(x[2])*(sin(ξ)*dξ[1]^3 - 3*cos(ξ)*dξ[1]*dξ[2] - sin(ξ)*dξ[3])) +
10*(-2*cos(x[2])*sin(ξ)*dx[2]*dξ[1] + cos(ξ)*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) + sin(x[2])*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2]))*(3*dx[1]*d2x[1]*dρ[2] + dρ[1]*d3x[1] +
dx[1]^3*dρ[3]) + 5*dx[1]*dρ[1]*(6*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2]) - 4*sin(ξ)*dξ[1]*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] +
cos(x[2])*d3x[2]) + 4*cos(x[2])*dx[2]*(sin(ξ)*dξ[1]^3 - 3*cos(ξ)*dξ[1]*dξ[2] - sin(ξ)*dξ[3]) + cos(ξ)*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 -
4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + sin(x[2])*(cos(ξ)*dξ[1]^4 + 6*sin(ξ)*dξ[1]^2*dξ[2] - 3*cos(ξ)*dξ[2]^2 - 4*cos(ξ)*dξ[1]*dξ[3] - sin(ξ)*dξ[4])) +
5*(cos(x[2])*cos(ξ)*dx[2] - sin(x[2])*sin(ξ)*dξ[1])*(3*d2x[1]^2*dρ[2] + 4*dx[1]*dρ[2]*d3x[1] + 6*dx[1]^2*d2x[1]*dρ[3] + dρ[1]*d4x[1] + dx[1]^4*dρ[4]) +
ρ*(10*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2])*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) + 10*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(sin(ξ)*dξ[1]^3 -
3*cos(ξ)*dξ[1]*dξ[2] - sin(ξ)*dξ[3]) - 5*sin(ξ)*dξ[1]*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) +
5*cos(x[2])*dx[2]*(cos(ξ)*dξ[1]^4 + 6*sin(ξ)*dξ[1]^2*dξ[2] - 3*cos(ξ)*dξ[2]^2 - 4*cos(ξ)*dξ[1]*dξ[3] - sin(ξ)*dξ[4]) + cos(ξ)*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] -
15*cos(x[2])*dx[2]*d2x[2]^2 - 10*cos(x[2])*dx[2]^2*d3x[2] - 10*sin(x[2])*d2x[2]*d3x[2] - 5*sin(x[2])*dx[2]*d4x[2] + cos(x[2])*d5x[2]) + sin(x[2])*(-(sin(ξ)*dξ[1]^5) +
10*cos(ξ)*dξ[1]^3*dξ[2] + 15*sin(ξ)*dξ[1]*dξ[2]^2 + 10*sin(ξ)*dξ[1]^2*dξ[3] - 10*cos(ξ)*dξ[2]*dξ[3] - 5*cos(ξ)*dξ[1]*dξ[4] - sin(ξ)*dξ[5])) +
cos(ξ)*sin(x[2])*(10*d2x[1]*dρ[2]*d3x[1] + 15*dx[1]*d2x[1]^2*dρ[3] + 10*dx[1]^2*d3x[1]*dρ[3] + 5*dx[1]*dρ[2]*d4x[1] + 10*dx[1]^3*d2x[1]*dρ[4] + dρ[1]*d5x[1] + dx[1]^5*dρ[5])

dxH6(x::AbstractVector{Float64}, dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, d4x::AbstractVector{Float64}, d5x::AbstractVector{Float64}, d6x::AbstractVector{Float64}, ξ::Float64, dξ::AbstractVector{Float64}, ρ::Float64, dρ::AbstractVector{Float64}) = 20*(-3*sin(ξ)*dξ[1]*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) + 3*cos(x[2])*dx[2]*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2]) +
cos(ξ)*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) + sin(x[2])*(sin(ξ)*dξ[1]^3 - 3*cos(ξ)*dξ[1]*dξ[2] - sin(ξ)*dξ[3]))*(3*dx[1]*d2x[1]*dρ[2] +
dρ[1]*d3x[1] + dx[1]^3*dρ[3]) + 15*(dρ[1]*d2x[1] + dx[1]^2*dρ[2])*(6*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2]) -
4*sin(ξ)*dξ[1]*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) + 4*cos(x[2])*dx[2]*(sin(ξ)*dξ[1]^3 - 3*cos(ξ)*dξ[1]*dξ[2] - sin(ξ)*dξ[3]) +
cos(ξ)*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + sin(x[2])*(cos(ξ)*dξ[1]^4 +
6*sin(ξ)*dξ[1]^2*dξ[2] - 3*cos(ξ)*dξ[2]^2 - 4*cos(ξ)*dξ[1]*dξ[3] - sin(ξ)*dξ[4])) + 15*(-2*cos(x[2])*sin(ξ)*dx[2]*dξ[1] + cos(ξ)*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) +
sin(x[2])*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2]))*(3*d2x[1]^2*dρ[2] + 4*dx[1]*dρ[2]*d3x[1] + 6*dx[1]^2*d2x[1]*dρ[3] + dρ[1]*d4x[1] + dx[1]^4*dρ[4]) +
6*dx[1]*dρ[1]*(10*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2])*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) + 10*(-(sin(x[2])*dx[2]^2) +
cos(x[2])*d2x[2])*(sin(ξ)*dξ[1]^3 - 3*cos(ξ)*dξ[1]*dξ[2] - sin(ξ)*dξ[3]) - 5*sin(ξ)*dξ[1]*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 -
4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + 5*cos(x[2])*dx[2]*(cos(ξ)*dξ[1]^4 + 6*sin(ξ)*dξ[1]^2*dξ[2] - 3*cos(ξ)*dξ[2]^2 - 4*cos(ξ)*dξ[1]*dξ[3] - sin(ξ)*dξ[4]) +
cos(ξ)*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] - 15*cos(x[2])*dx[2]*d2x[2]^2 - 10*cos(x[2])*dx[2]^2*d3x[2] - 10*sin(x[2])*d2x[2]*d3x[2] - 5*sin(x[2])*dx[2]*d4x[2] +
cos(x[2])*d5x[2]) + sin(x[2])*(-(sin(ξ)*dξ[1]^5) + 10*cos(ξ)*dξ[1]^3*dξ[2] + 15*sin(ξ)*dξ[1]*dξ[2]^2 + 10*sin(ξ)*dξ[1]^2*dξ[3] - 10*cos(ξ)*dξ[2]*dξ[3] - 5*cos(ξ)*dξ[1]*dξ[4] -
sin(ξ)*dξ[5])) + 6*(cos(x[2])*cos(ξ)*dx[2] - sin(x[2])*sin(ξ)*dξ[1])*(10*d2x[1]*dρ[2]*d3x[1] + 15*dx[1]*d2x[1]^2*dρ[3] + 10*dx[1]^2*d3x[1]*dρ[3] + 5*dx[1]*dρ[2]*d4x[1] +
10*dx[1]^3*d2x[1]*dρ[4] + dρ[1]*d5x[1] + dx[1]^5*dρ[5]) + ρ*(20*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2])*(sin(ξ)*dξ[1]^3 - 3*cos(ξ)*dξ[1]*dξ[2] -
sin(ξ)*dξ[3]) + 15*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2])*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) +
15*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(cos(ξ)*dξ[1]^4 + 6*sin(ξ)*dξ[1]^2*dξ[2] - 3*cos(ξ)*dξ[2]^2 - 4*cos(ξ)*dξ[1]*dξ[3] - sin(ξ)*dξ[4]) - 6*sin(ξ)*dξ[1]*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] - 15*cos(x[2])*dx[2]*d2x[2]^2 - 10*cos(x[2])*dx[2]^2*d3x[2] - 10*sin(x[2])*d2x[2]*d3x[2] - 5*sin(x[2])*dx[2]*d4x[2] + cos(x[2])*d5x[2]) +
6*cos(x[2])*dx[2]*(-(sin(ξ)*dξ[1]^5) + 10*cos(ξ)*dξ[1]^3*dξ[2] + 15*sin(ξ)*dξ[1]*dξ[2]^2 + 10*sin(ξ)*dξ[1]^2*dξ[3] - 10*cos(ξ)*dξ[2]*dξ[3] - 5*cos(ξ)*dξ[1]*dξ[4] - sin(ξ)*dξ[5]) +
cos(ξ)*(-(sin(x[2])*dx[2]^6) + 15*cos(x[2])*dx[2]^4*d2x[2] + 45*sin(x[2])*dx[2]^2*d2x[2]^2 - 15*cos(x[2])*d2x[2]^3 + 20*sin(x[2])*dx[2]^3*d3x[2] -
60*cos(x[2])*dx[2]*d2x[2]*d3x[2] - 10*sin(x[2])*d3x[2]^2 - 15*cos(x[2])*dx[2]^2*d4x[2] - 15*sin(x[2])*d2x[2]*d4x[2] - 6*sin(x[2])*dx[2]*d5x[2] + cos(x[2])*d6x[2]) +
sin(x[2])*(-(cos(ξ)*dξ[1]^6) - 15*sin(ξ)*dξ[1]^4*dξ[2] + 45*cos(ξ)*dξ[1]^2*dξ[2]^2 + 15*sin(ξ)*dξ[2]^3 + 20*cos(ξ)*dξ[1]^3*dξ[3] + 60*sin(ξ)*dξ[1]*dξ[2]*dξ[3] -
10*cos(ξ)*dξ[3]^2 + 15*sin(ξ)*dξ[1]^2*dξ[4] - 15*cos(ξ)*dξ[2]*dξ[4] - 6*cos(ξ)*dξ[1]*dξ[5] - sin(ξ)*dξ[6])) + cos(ξ)*sin(x[2])*(10*dρ[2]*d3x[1]^2 + 15*d2x[1]^3*dρ[3] +
60*dx[1]*d2x[1]*d3x[1]*dρ[3] + 15*d2x[1]*dρ[2]*d4x[1] + 15*dx[1]^2*dρ[3]*d4x[1] + 45*dx[1]^2*d2x[1]^2*dρ[4] + 20*dx[1]^3*d3x[1]*dρ[4] + 6*dx[1]*dρ[2]*d5x[1] +
15*dx[1]^4*d2x[1]*dρ[5] + dρ[1]*d6x[1] + dx[1]^6*dρ[6])

dxH7(x::AbstractVector{Float64}, dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, d4x::AbstractVector{Float64}, d5x::AbstractVector{Float64}, d6x::AbstractVector{Float64}, d7x::AbstractVector{Float64}, ξ::Float64, dξ::AbstractVector{Float64}, ρ::Float64, dρ::AbstractVector{Float64}) = 35*(3*dx[1]*d2x[1]*dρ[2] + dρ[1]*d3x[1] + dx[1]^3*dρ[3])*(6*(-(sin(x[2])*dx[2]^2) +
cos(x[2])*d2x[2])*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2]) - 4*sin(ξ)*dξ[1]*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) +
4*cos(x[2])*dx[2]*(sin(ξ)*dξ[1]^3 - 3*cos(ξ)*dξ[1]*dξ[2] - sin(ξ)*dξ[3]) + cos(ξ)*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 -
4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + sin(x[2])*(cos(ξ)*dξ[1]^4 + 6*sin(ξ)*dξ[1]^2*dξ[2] - 3*cos(ξ)*dξ[2]^2 - 4*cos(ξ)*dξ[1]*dξ[3] - sin(ξ)*dξ[4])) +
35*(-3*sin(ξ)*dξ[1]*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) + 3*cos(x[2])*dx[2]*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2]) + cos(ξ)*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] +
cos(x[2])*d3x[2]) + sin(x[2])*(sin(ξ)*dξ[1]^3 - 3*cos(ξ)*dξ[1]*dξ[2] - sin(ξ)*dξ[3]))*(3*d2x[1]^2*dρ[2] + 4*dx[1]*dρ[2]*d3x[1] + 6*dx[1]^2*d2x[1]*dρ[3] + dρ[1]*d4x[1] +
dx[1]^4*dρ[4]) + 21*(dρ[1]*d2x[1] + dx[1]^2*dρ[2])*(10*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2])*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) +
10*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(sin(ξ)*dξ[1]^3 - 3*cos(ξ)*dξ[1]*dξ[2] - sin(ξ)*dξ[3]) - 5*sin(ξ)*dξ[1]*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] -
3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + 5*cos(x[2])*dx[2]*(cos(ξ)*dξ[1]^4 + 6*sin(ξ)*dξ[1]^2*dξ[2] - 3*cos(ξ)*dξ[2]^2 - 4*cos(ξ)*dξ[1]*dξ[3] -
sin(ξ)*dξ[4]) + cos(ξ)*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] - 15*cos(x[2])*dx[2]*d2x[2]^2 - 10*cos(x[2])*dx[2]^2*d3x[2] - 10*sin(x[2])*d2x[2]*d3x[2] -
5*sin(x[2])*dx[2]*d4x[2] + cos(x[2])*d5x[2]) + sin(x[2])*(-(sin(ξ)*dξ[1]^5) + 10*cos(ξ)*dξ[1]^3*dξ[2] + 15*sin(ξ)*dξ[1]*dξ[2]^2 + 10*sin(ξ)*dξ[1]^2*dξ[3] -
10*cos(ξ)*dξ[2]*dξ[3] - 5*cos(ξ)*dξ[1]*dξ[4] - sin(ξ)*dξ[5])) + 21*(-2*cos(x[2])*sin(ξ)*dx[2]*dξ[1] + cos(ξ)*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) + sin(x[2])*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2]))*(10*d2x[1]*dρ[2]*d3x[1] + 15*dx[1]*d2x[1]^2*dρ[3] + 10*dx[1]^2*d3x[1]*dρ[3] + 5*dx[1]*dρ[2]*d4x[1] + 10*dx[1]^3*d2x[1]*dρ[4] + dρ[1]*d5x[1] +
dx[1]^5*dρ[5]) + 7*dx[1]*dρ[1]*(20*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2])*(sin(ξ)*dξ[1]^3 - 3*cos(ξ)*dξ[1]*dξ[2] - sin(ξ)*dξ[3]) + 15*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2])*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + 15*(-(sin(x[2])*dx[2]^2) +
cos(x[2])*d2x[2])*(cos(ξ)*dξ[1]^4 + 6*sin(ξ)*dξ[1]^2*dξ[2] - 3*cos(ξ)*dξ[2]^2 - 4*cos(ξ)*dξ[1]*dξ[3] - sin(ξ)*dξ[4]) - 6*sin(ξ)*dξ[1]*(cos(x[2])*dx[2]^5 +
10*sin(x[2])*dx[2]^3*d2x[2] - 15*cos(x[2])*dx[2]*d2x[2]^2 - 10*cos(x[2])*dx[2]^2*d3x[2] - 10*sin(x[2])*d2x[2]*d3x[2] - 5*sin(x[2])*dx[2]*d4x[2] + cos(x[2])*d5x[2]) +
6*cos(x[2])*dx[2]*(-(sin(ξ)*dξ[1]^5) + 10*cos(ξ)*dξ[1]^3*dξ[2] + 15*sin(ξ)*dξ[1]*dξ[2]^2 + 10*sin(ξ)*dξ[1]^2*dξ[3] - 10*cos(ξ)*dξ[2]*dξ[3] - 5*cos(ξ)*dξ[1]*dξ[4] -
sin(ξ)*dξ[5]) + cos(ξ)*(-(sin(x[2])*dx[2]^6) + 15*cos(x[2])*dx[2]^4*d2x[2] + 45*sin(x[2])*dx[2]^2*d2x[2]^2 - 15*cos(x[2])*d2x[2]^3 + 20*sin(x[2])*dx[2]^3*d3x[2] -
60*cos(x[2])*dx[2]*d2x[2]*d3x[2] - 10*sin(x[2])*d3x[2]^2 - 15*cos(x[2])*dx[2]^2*d4x[2] - 15*sin(x[2])*d2x[2]*d4x[2] - 6*sin(x[2])*dx[2]*d5x[2] + cos(x[2])*d6x[2]) +
sin(x[2])*(-(cos(ξ)*dξ[1]^6) - 15*sin(ξ)*dξ[1]^4*dξ[2] + 45*cos(ξ)*dξ[1]^2*dξ[2]^2 + 15*sin(ξ)*dξ[2]^3 + 20*cos(ξ)*dξ[1]^3*dξ[3] + 60*sin(ξ)*dξ[1]*dξ[2]*dξ[3] -
10*cos(ξ)*dξ[3]^2 + 15*sin(ξ)*dξ[1]^2*dξ[4] - 15*cos(ξ)*dξ[2]*dξ[4] - 6*cos(ξ)*dξ[1]*dξ[5] - sin(ξ)*dξ[6])) + 7*(cos(x[2])*cos(ξ)*dx[2] -
sin(x[2])*sin(ξ)*dξ[1])*(10*dρ[2]*d3x[1]^2 + 15*d2x[1]^3*dρ[3] + 60*dx[1]*d2x[1]*d3x[1]*dρ[3] + 15*d2x[1]*dρ[2]*d4x[1] + 15*dx[1]^2*dρ[3]*d4x[1] +
45*dx[1]^2*d2x[1]^2*dρ[4] + 20*dx[1]^3*d3x[1]*dρ[4] + 6*dx[1]*dρ[2]*d5x[1] + 15*dx[1]^4*d2x[1]*dρ[5] + dρ[1]*d6x[1] + dx[1]^6*dρ[6]) + ρ*(35*(sin(ξ)*dξ[1]^3 -
3*cos(ξ)*dξ[1]*dξ[2] - sin(ξ)*dξ[3])*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) +
35*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2])*(cos(ξ)*dξ[1]^4 + 6*sin(ξ)*dξ[1]^2*dξ[2] - 3*cos(ξ)*dξ[2]^2 - 4*cos(ξ)*dξ[1]*dξ[3] -
sin(ξ)*dξ[4]) + 21*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2])*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] - 15*cos(x[2])*dx[2]*d2x[2]^2 - 10*cos(x[2])*dx[2]^2*d3x[2] -
10*sin(x[2])*d2x[2]*d3x[2] - 5*sin(x[2])*dx[2]*d4x[2] + cos(x[2])*d5x[2]) + 21*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(-(sin(ξ)*dξ[1]^5) + 10*cos(ξ)*dξ[1]^3*dξ[2] +
15*sin(ξ)*dξ[1]*dξ[2]^2 + 10*sin(ξ)*dξ[1]^2*dξ[3] - 10*cos(ξ)*dξ[2]*dξ[3] - 5*cos(ξ)*dξ[1]*dξ[4] - sin(ξ)*dξ[5]) - 7*sin(ξ)*dξ[1]*(-(sin(x[2])*dx[2]^6) +
15*cos(x[2])*dx[2]^4*d2x[2] + 45*sin(x[2])*dx[2]^2*d2x[2]^2 - 15*cos(x[2])*d2x[2]^3 + 20*sin(x[2])*dx[2]^3*d3x[2] - 60*cos(x[2])*dx[2]*d2x[2]*d3x[2] - 10*sin(x[2])*d3x[2]^2 -
15*cos(x[2])*dx[2]^2*d4x[2] - 15*sin(x[2])*d2x[2]*d4x[2] - 6*sin(x[2])*dx[2]*d5x[2] + cos(x[2])*d6x[2]) + 7*cos(x[2])*dx[2]*(-(cos(ξ)*dξ[1]^6) - 15*sin(ξ)*dξ[1]^4*dξ[2] +
45*cos(ξ)*dξ[1]^2*dξ[2]^2 + 15*sin(ξ)*dξ[2]^3 + 20*cos(ξ)*dξ[1]^3*dξ[3] + 60*sin(ξ)*dξ[1]*dξ[2]*dξ[3] - 10*cos(ξ)*dξ[3]^2 + 15*sin(ξ)*dξ[1]^2*dξ[4] - 15*cos(ξ)*dξ[2]*dξ[4] -
6*cos(ξ)*dξ[1]*dξ[5] - sin(ξ)*dξ[6]) + cos(ξ)*(-(cos(x[2])*dx[2]^7) - 21*sin(x[2])*dx[2]^5*d2x[2] + 105*cos(x[2])*dx[2]^3*d2x[2]^2 + 105*sin(x[2])*dx[2]*d2x[2]^3 +
35*cos(x[2])*dx[2]^4*d3x[2] + 210*sin(x[2])*dx[2]^2*d2x[2]*d3x[2] - 105*cos(x[2])*d2x[2]^2*d3x[2] - 70*cos(x[2])*dx[2]*d3x[2]^2 + 35*sin(x[2])*dx[2]^3*d4x[2] -
105*cos(x[2])*dx[2]*d2x[2]*d4x[2] - 35*sin(x[2])*d3x[2]*d4x[2] - 21*cos(x[2])*dx[2]^2*d5x[2] - 21*sin(x[2])*d2x[2]*d5x[2] - 7*sin(x[2])*dx[2]*d6x[2] + cos(x[2])*d7x[2]) +
sin(x[2])*(sin(ξ)*dξ[1]^7 - 21*cos(ξ)*dξ[1]^5*dξ[2] - 105*sin(ξ)*dξ[1]^3*dξ[2]^2 + 105*cos(ξ)*dξ[1]*dξ[2]^3 - 35*sin(ξ)*dξ[1]^4*dξ[3] + 210*cos(ξ)*dξ[1]^2*dξ[2]*dξ[3] +
105*sin(ξ)*dξ[2]^2*dξ[3] + 70*sin(ξ)*dξ[1]*dξ[3]^2 + 35*cos(ξ)*dξ[1]^3*dξ[4] + 105*sin(ξ)*dξ[1]*dξ[2]*dξ[4] - 35*cos(ξ)*dξ[3]*dξ[4] + 21*sin(ξ)*dξ[1]^2*dξ[5] -
21*cos(ξ)*dξ[2]*dξ[5] - 7*cos(ξ)*dξ[1]*dξ[6] - sin(ξ)*dξ[7])) + cos(ξ)*sin(x[2])*(105*d2x[1]^2*d3x[1]*dρ[3] + 70*dx[1]*d3x[1]^2*dρ[3] + 35*dρ[2]*d3x[1]*d4x[1] +
105*dx[1]*d2x[1]*dρ[3]*d4x[1] + 105*dx[1]*d2x[1]^3*dρ[4] + 210*dx[1]^2*d2x[1]*d3x[1]*dρ[4] + 35*dx[1]^3*d4x[1]*dρ[4] + 21*d2x[1]*dρ[2]*d5x[1] + 21*dx[1]^2*dρ[3]*d5x[1] +
105*dx[1]^3*d2x[1]^2*dρ[5] + 35*dx[1]^4*d3x[1]*dρ[5] + 7*dx[1]*dρ[2]*d6x[1] + 21*dx[1]^5*d2x[1]*dρ[6] + dρ[1]*d7x[1] + dx[1]^7*dρ[7])

dxH8(x::AbstractVector{Float64}, dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, d4x::AbstractVector{Float64}, d5x::AbstractVector{Float64}, d6x::AbstractVector{Float64}, d7x::AbstractVector{Float64}, d8x::AbstractVector{Float64}, ξ::Float64, dξ::AbstractVector{Float64}, ρ::Float64, dρ::AbstractVector{Float64}) = 70*(6*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2]) -
4*sin(ξ)*dξ[1]*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) + 4*cos(x[2])*dx[2]*(sin(ξ)*dξ[1]^3 - 3*cos(ξ)*dξ[1]*dξ[2] - sin(ξ)*dξ[3]) +
cos(ξ)*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + sin(x[2])*(cos(ξ)*dξ[1]^4 +
6*sin(ξ)*dξ[1]^2*dξ[2] - 3*cos(ξ)*dξ[2]^2 - 4*cos(ξ)*dξ[1]*dξ[3] - sin(ξ)*dξ[4]))*(3*d2x[1]^2*dρ[2] + 4*dx[1]*dρ[2]*d3x[1] + 6*dx[1]^2*d2x[1]*dρ[3] + dρ[1]*d4x[1] +
dx[1]^4*dρ[4]) + 56*(3*dx[1]*d2x[1]*dρ[2] + dρ[1]*d3x[1] + dx[1]^3*dρ[3])*(10*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2])*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] +
cos(x[2])*d3x[2]) + 10*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(sin(ξ)*dξ[1]^3 - 3*cos(ξ)*dξ[1]*dξ[2] - sin(ξ)*dξ[3]) - 5*sin(ξ)*dξ[1]*(sin(x[2])*dx[2]^4 -
6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + 5*cos(x[2])*dx[2]*(cos(ξ)*dξ[1]^4 + 6*sin(ξ)*dξ[1]^2*dξ[2] -
3*cos(ξ)*dξ[2]^2 - 4*cos(ξ)*dξ[1]*dξ[3] - sin(ξ)*dξ[4]) + cos(ξ)*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] - 15*cos(x[2])*dx[2]*d2x[2]^2 -
10*cos(x[2])*dx[2]^2*d3x[2] - 10*sin(x[2])*d2x[2]*d3x[2] - 5*sin(x[2])*dx[2]*d4x[2] + cos(x[2])*d5x[2]) + sin(x[2])*(-(sin(ξ)*dξ[1]^5) + 10*cos(ξ)*dξ[1]^3*dξ[2] +
15*sin(ξ)*dξ[1]*dξ[2]^2 + 10*sin(ξ)*dξ[1]^2*dξ[3] - 10*cos(ξ)*dξ[2]*dξ[3] - 5*cos(ξ)*dξ[1]*dξ[4] - sin(ξ)*dξ[5])) + 56*(-3*sin(ξ)*dξ[1]*(-(sin(x[2])*dx[2]^2) +
cos(x[2])*d2x[2]) + 3*cos(x[2])*dx[2]*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2]) + cos(ξ)*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) +
sin(x[2])*(sin(ξ)*dξ[1]^3 - 3*cos(ξ)*dξ[1]*dξ[2] - sin(ξ)*dξ[3]))*(10*d2x[1]*dρ[2]*d3x[1] + 15*dx[1]*d2x[1]^2*dρ[3] + 10*dx[1]^2*d3x[1]*dρ[3] + 5*dx[1]*dρ[2]*d4x[1] +
10*dx[1]^3*d2x[1]*dρ[4] + dρ[1]*d5x[1] + dx[1]^5*dρ[5]) + 28*(dρ[1]*d2x[1] + dx[1]^2*dρ[2])*(20*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] +
cos(x[2])*d3x[2])*(sin(ξ)*dξ[1]^3 - 3*cos(ξ)*dξ[1]*dξ[2] - sin(ξ)*dξ[3]) + 15*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2])*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] -
3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + 15*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(cos(ξ)*dξ[1]^4 + 6*sin(ξ)*dξ[1]^2*dξ[2] -
3*cos(ξ)*dξ[2]^2 - 4*cos(ξ)*dξ[1]*dξ[3] - sin(ξ)*dξ[4]) - 6*sin(ξ)*dξ[1]*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] - 15*cos(x[2])*dx[2]*d2x[2]^2 -
10*cos(x[2])*dx[2]^2*d3x[2] - 10*sin(x[2])*d2x[2]*d3x[2] - 5*sin(x[2])*dx[2]*d4x[2] + cos(x[2])*d5x[2]) + 6*cos(x[2])*dx[2]*(-(sin(ξ)*dξ[1]^5) + 10*cos(ξ)*dξ[1]^3*dξ[2] +
15*sin(ξ)*dξ[1]*dξ[2]^2 + 10*sin(ξ)*dξ[1]^2*dξ[3] - 10*cos(ξ)*dξ[2]*dξ[3] - 5*cos(ξ)*dξ[1]*dξ[4] - sin(ξ)*dξ[5]) + cos(ξ)*(-(sin(x[2])*dx[2]^6) +
15*cos(x[2])*dx[2]^4*d2x[2] + 45*sin(x[2])*dx[2]^2*d2x[2]^2 - 15*cos(x[2])*d2x[2]^3 + 20*sin(x[2])*dx[2]^3*d3x[2] - 60*cos(x[2])*dx[2]*d2x[2]*d3x[2] - 10*sin(x[2])*d3x[2]^2 -
15*cos(x[2])*dx[2]^2*d4x[2] - 15*sin(x[2])*d2x[2]*d4x[2] - 6*sin(x[2])*dx[2]*d5x[2] + cos(x[2])*d6x[2]) + sin(x[2])*(-(cos(ξ)*dξ[1]^6) - 15*sin(ξ)*dξ[1]^4*dξ[2] +
45*cos(ξ)*dξ[1]^2*dξ[2]^2 + 15*sin(ξ)*dξ[2]^3 + 20*cos(ξ)*dξ[1]^3*dξ[3] + 60*sin(ξ)*dξ[1]*dξ[2]*dξ[3] - 10*cos(ξ)*dξ[3]^2 + 15*sin(ξ)*dξ[1]^2*dξ[4] - 15*cos(ξ)*dξ[2]*dξ[4] -
6*cos(ξ)*dξ[1]*dξ[5] - sin(ξ)*dξ[6])) + 28*(-2*cos(x[2])*sin(ξ)*dx[2]*dξ[1] + cos(ξ)*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) + sin(x[2])*(-(cos(ξ)*dξ[1]^2) -
sin(ξ)*dξ[2]))*(10*dρ[2]*d3x[1]^2 + 15*d2x[1]^3*dρ[3] + 60*dx[1]*d2x[1]*d3x[1]*dρ[3] + 15*d2x[1]*dρ[2]*d4x[1] + 15*dx[1]^2*dρ[3]*d4x[1] + 45*dx[1]^2*d2x[1]^2*dρ[4] +
20*dx[1]^3*d3x[1]*dρ[4] + 6*dx[1]*dρ[2]*d5x[1] + 15*dx[1]^4*d2x[1]*dρ[5] + dρ[1]*d6x[1] + dx[1]^6*dρ[6]) + 8*dx[1]*dρ[1]*(35*(sin(ξ)*dξ[1]^3 - 3*cos(ξ)*dξ[1]*dξ[2] -
sin(ξ)*dξ[3])*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + 35*(-(cos(x[2])*dx[2]^3) -
3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2])*(cos(ξ)*dξ[1]^4 + 6*sin(ξ)*dξ[1]^2*dξ[2] - 3*cos(ξ)*dξ[2]^2 - 4*cos(ξ)*dξ[1]*dξ[3] - sin(ξ)*dξ[4]) + 21*(-(cos(ξ)*dξ[1]^2) -
sin(ξ)*dξ[2])*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] - 15*cos(x[2])*dx[2]*d2x[2]^2 - 10*cos(x[2])*dx[2]^2*d3x[2] - 10*sin(x[2])*d2x[2]*d3x[2] -
5*sin(x[2])*dx[2]*d4x[2] + cos(x[2])*d5x[2]) + 21*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(-(sin(ξ)*dξ[1]^5) + 10*cos(ξ)*dξ[1]^3*dξ[2] + 15*sin(ξ)*dξ[1]*dξ[2]^2 +
10*sin(ξ)*dξ[1]^2*dξ[3] - 10*cos(ξ)*dξ[2]*dξ[3] - 5*cos(ξ)*dξ[1]*dξ[4] - sin(ξ)*dξ[5]) - 7*sin(ξ)*dξ[1]*(-(sin(x[2])*dx[2]^6) + 15*cos(x[2])*dx[2]^4*d2x[2] +
45*sin(x[2])*dx[2]^2*d2x[2]^2 - 15*cos(x[2])*d2x[2]^3 + 20*sin(x[2])*dx[2]^3*d3x[2] - 60*cos(x[2])*dx[2]*d2x[2]*d3x[2] - 10*sin(x[2])*d3x[2]^2 - 15*cos(x[2])*dx[2]^2*d4x[2] -
15*sin(x[2])*d2x[2]*d4x[2] - 6*sin(x[2])*dx[2]*d5x[2] + cos(x[2])*d6x[2]) + 7*cos(x[2])*dx[2]*(-(cos(ξ)*dξ[1]^6) - 15*sin(ξ)*dξ[1]^4*dξ[2] + 45*cos(ξ)*dξ[1]^2*dξ[2]^2 +
15*sin(ξ)*dξ[2]^3 + 20*cos(ξ)*dξ[1]^3*dξ[3] + 60*sin(ξ)*dξ[1]*dξ[2]*dξ[3] - 10*cos(ξ)*dξ[3]^2 + 15*sin(ξ)*dξ[1]^2*dξ[4] - 15*cos(ξ)*dξ[2]*dξ[4] - 6*cos(ξ)*dξ[1]*dξ[5] -
sin(ξ)*dξ[6]) + cos(ξ)*(-(cos(x[2])*dx[2]^7) - 21*sin(x[2])*dx[2]^5*d2x[2] + 105*cos(x[2])*dx[2]^3*d2x[2]^2 + 105*sin(x[2])*dx[2]*d2x[2]^3 + 35*cos(x[2])*dx[2]^4*d3x[2] +
210*sin(x[2])*dx[2]^2*d2x[2]*d3x[2] - 105*cos(x[2])*d2x[2]^2*d3x[2] - 70*cos(x[2])*dx[2]*d3x[2]^2 + 35*sin(x[2])*dx[2]^3*d4x[2] - 105*cos(x[2])*dx[2]*d2x[2]*d4x[2] -
35*sin(x[2])*d3x[2]*d4x[2] - 21*cos(x[2])*dx[2]^2*d5x[2] - 21*sin(x[2])*d2x[2]*d5x[2] - 7*sin(x[2])*dx[2]*d6x[2] + cos(x[2])*d7x[2]) + sin(x[2])*(sin(ξ)*dξ[1]^7 -
21*cos(ξ)*dξ[1]^5*dξ[2] - 105*sin(ξ)*dξ[1]^3*dξ[2]^2 + 105*cos(ξ)*dξ[1]*dξ[2]^3 - 35*sin(ξ)*dξ[1]^4*dξ[3] + 210*cos(ξ)*dξ[1]^2*dξ[2]*dξ[3] + 105*sin(ξ)*dξ[2]^2*dξ[3] +
70*sin(ξ)*dξ[1]*dξ[3]^2 + 35*cos(ξ)*dξ[1]^3*dξ[4] + 105*sin(ξ)*dξ[1]*dξ[2]*dξ[4] - 35*cos(ξ)*dξ[3]*dξ[4] + 21*sin(ξ)*dξ[1]^2*dξ[5] - 21*cos(ξ)*dξ[2]*dξ[5] -
7*cos(ξ)*dξ[1]*dξ[6] - sin(ξ)*dξ[7])) + 8*(cos(x[2])*cos(ξ)*dx[2] - sin(x[2])*sin(ξ)*dξ[1])*(105*d2x[1]^2*d3x[1]*dρ[3] + 70*dx[1]*d3x[1]^2*dρ[3] + 35*dρ[2]*d3x[1]*d4x[1] +
105*dx[1]*d2x[1]*dρ[3]*d4x[1] + 105*dx[1]*d2x[1]^3*dρ[4] + 210*dx[1]^2*d2x[1]*d3x[1]*dρ[4] + 35*dx[1]^3*d4x[1]*dρ[4] + 21*d2x[1]*dρ[2]*d5x[1] + 21*dx[1]^2*dρ[3]*d5x[1] +
105*dx[1]^3*d2x[1]^2*dρ[5] + 35*dx[1]^4*d3x[1]*dρ[5] + 7*dx[1]*dρ[2]*d6x[1] + 21*dx[1]^5*d2x[1]*dρ[6] + dρ[1]*d7x[1] + dx[1]^7*dρ[7]) + ρ*(70*(sin(x[2])*dx[2]^4 -
6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2])*(cos(ξ)*dξ[1]^4 + 6*sin(ξ)*dξ[1]^2*dξ[2] - 3*cos(ξ)*dξ[2]^2 -
4*cos(ξ)*dξ[1]*dξ[3] - sin(ξ)*dξ[4]) + 56*(sin(ξ)*dξ[1]^3 - 3*cos(ξ)*dξ[1]*dξ[2] - sin(ξ)*dξ[3])*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] -
15*cos(x[2])*dx[2]*d2x[2]^2 - 10*cos(x[2])*dx[2]^2*d3x[2] - 10*sin(x[2])*d2x[2]*d3x[2] - 5*sin(x[2])*dx[2]*d4x[2] + cos(x[2])*d5x[2]) + 56*(-(cos(x[2])*dx[2]^3) -
3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2])*(-(sin(ξ)*dξ[1]^5) + 10*cos(ξ)*dξ[1]^3*dξ[2] + 15*sin(ξ)*dξ[1]*dξ[2]^2 + 10*sin(ξ)*dξ[1]^2*dξ[3] - 10*cos(ξ)*dξ[2]*dξ[3] -
5*cos(ξ)*dξ[1]*dξ[4] - sin(ξ)*dξ[5]) + 28*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2])*(-(sin(x[2])*dx[2]^6) + 15*cos(x[2])*dx[2]^4*d2x[2] + 45*sin(x[2])*dx[2]^2*d2x[2]^2 -
15*cos(x[2])*d2x[2]^3 + 20*sin(x[2])*dx[2]^3*d3x[2] - 60*cos(x[2])*dx[2]*d2x[2]*d3x[2] - 10*sin(x[2])*d3x[2]^2 - 15*cos(x[2])*dx[2]^2*d4x[2] - 15*sin(x[2])*d2x[2]*d4x[2] -
6*sin(x[2])*dx[2]*d5x[2] + cos(x[2])*d6x[2]) + 28*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(-(cos(ξ)*dξ[1]^6) - 15*sin(ξ)*dξ[1]^4*dξ[2] + 45*cos(ξ)*dξ[1]^2*dξ[2]^2 +
15*sin(ξ)*dξ[2]^3 + 20*cos(ξ)*dξ[1]^3*dξ[3] + 60*sin(ξ)*dξ[1]*dξ[2]*dξ[3] - 10*cos(ξ)*dξ[3]^2 + 15*sin(ξ)*dξ[1]^2*dξ[4] - 15*cos(ξ)*dξ[2]*dξ[4] - 6*cos(ξ)*dξ[1]*dξ[5] -
sin(ξ)*dξ[6]) - 8*sin(ξ)*dξ[1]*(-(cos(x[2])*dx[2]^7) - 21*sin(x[2])*dx[2]^5*d2x[2] + 105*cos(x[2])*dx[2]^3*d2x[2]^2 + 105*sin(x[2])*dx[2]*d2x[2]^3 +
35*cos(x[2])*dx[2]^4*d3x[2] + 210*sin(x[2])*dx[2]^2*d2x[2]*d3x[2] - 105*cos(x[2])*d2x[2]^2*d3x[2] - 70*cos(x[2])*dx[2]*d3x[2]^2 + 35*sin(x[2])*dx[2]^3*d4x[2] -
105*cos(x[2])*dx[2]*d2x[2]*d4x[2] - 35*sin(x[2])*d3x[2]*d4x[2] - 21*cos(x[2])*dx[2]^2*d5x[2] - 21*sin(x[2])*d2x[2]*d5x[2] - 7*sin(x[2])*dx[2]*d6x[2] + cos(x[2])*d7x[2]) +
8*cos(x[2])*dx[2]*(sin(ξ)*dξ[1]^7 - 21*cos(ξ)*dξ[1]^5*dξ[2] - 105*sin(ξ)*dξ[1]^3*dξ[2]^2 + 105*cos(ξ)*dξ[1]*dξ[2]^3 - 35*sin(ξ)*dξ[1]^4*dξ[3] +
210*cos(ξ)*dξ[1]^2*dξ[2]*dξ[3] + 105*sin(ξ)*dξ[2]^2*dξ[3] + 70*sin(ξ)*dξ[1]*dξ[3]^2 + 35*cos(ξ)*dξ[1]^3*dξ[4] + 105*sin(ξ)*dξ[1]*dξ[2]*dξ[4] - 35*cos(ξ)*dξ[3]*dξ[4] +
21*sin(ξ)*dξ[1]^2*dξ[5] - 21*cos(ξ)*dξ[2]*dξ[5] - 7*cos(ξ)*dξ[1]*dξ[6] - sin(ξ)*dξ[7]) + cos(ξ)*(sin(x[2])*dx[2]^8 - 28*cos(x[2])*dx[2]^6*d2x[2] -
210*sin(x[2])*dx[2]^4*d2x[2]^2 + 420*cos(x[2])*dx[2]^2*d2x[2]^3 + 105*sin(x[2])*d2x[2]^4 - 56*sin(x[2])*dx[2]^5*d3x[2] + 560*cos(x[2])*dx[2]^3*d2x[2]*d3x[2] +
840*sin(x[2])*dx[2]*d2x[2]^2*d3x[2] + 280*sin(x[2])*dx[2]^2*d3x[2]^2 - 280*cos(x[2])*d2x[2]*d3x[2]^2 + 70*cos(x[2])*dx[2]^4*d4x[2] + 420*sin(x[2])*dx[2]^2*d2x[2]*d4x[2] -
210*cos(x[2])*d2x[2]^2*d4x[2] - 280*cos(x[2])*dx[2]*d3x[2]*d4x[2] - 35*sin(x[2])*d4x[2]^2 + 56*sin(x[2])*dx[2]^3*d5x[2] - 168*cos(x[2])*dx[2]*d2x[2]*d5x[2] -
56*sin(x[2])*d3x[2]*d5x[2] - 28*cos(x[2])*dx[2]^2*d6x[2] - 28*sin(x[2])*d2x[2]*d6x[2] - 8*sin(x[2])*dx[2]*d7x[2] + cos(x[2])*d8x[2]) + sin(x[2])*(cos(ξ)*dξ[1]^8 +
28*sin(ξ)*dξ[1]^6*dξ[2] - 210*cos(ξ)*dξ[1]^4*dξ[2]^2 - 420*sin(ξ)*dξ[1]^2*dξ[2]^3 + 105*cos(ξ)*dξ[2]^4 - 56*cos(ξ)*dξ[1]^5*dξ[3] - 560*sin(ξ)*dξ[1]^3*dξ[2]*dξ[3] +
840*cos(ξ)*dξ[1]*dξ[2]^2*dξ[3] + 280*cos(ξ)*dξ[1]^2*dξ[3]^2 + 280*sin(ξ)*dξ[2]*dξ[3]^2 - 70*sin(ξ)*dξ[1]^4*dξ[4] + 420*cos(ξ)*dξ[1]^2*dξ[2]*dξ[4] +
210*sin(ξ)*dξ[2]^2*dξ[4] + 280*sin(ξ)*dξ[1]*dξ[3]*dξ[4] - 35*cos(ξ)*dξ[4]^2 + 56*cos(ξ)*dξ[1]^3*dξ[5] + 168*sin(ξ)*dξ[1]*dξ[2]*dξ[5] - 56*cos(ξ)*dξ[3]*dξ[5] +
28*sin(ξ)*dξ[1]^2*dξ[6] - 28*cos(ξ)*dξ[2]*dξ[6] - 8*cos(ξ)*dξ[1]*dξ[7] - sin(ξ)*dξ[8])) + cos(ξ)*sin(x[2])*(280*d2x[1]*d3x[1]^2*dρ[3] + 210*d2x[1]^2*dρ[3]*d4x[1] +
280*dx[1]*d3x[1]*dρ[3]*d4x[1] + 35*dρ[2]*d4x[1]^2 + 105*d2x[1]^4*dρ[4] + 840*dx[1]*d2x[1]^2*d3x[1]*dρ[4] + 280*dx[1]^2*d3x[1]^2*dρ[4] + 420*dx[1]^2*d2x[1]*d4x[1]*dρ[4] +
56*dρ[2]*d3x[1]*d5x[1] + 168*dx[1]*d2x[1]*dρ[3]*d5x[1] + 56*dx[1]^3*dρ[4]*d5x[1] + 420*dx[1]^2*d2x[1]^3*dρ[5] + 560*dx[1]^3*d2x[1]*d3x[1]*dρ[5] + 70*dx[1]^4*d4x[1]*dρ[5] +
28*d2x[1]*dρ[2]*d6x[1] + 28*dx[1]^2*dρ[3]*d6x[1] + 210*dx[1]^4*d2x[1]^2*dρ[6] + 56*dx[1]^5*d3x[1]*dρ[6] + 8*dx[1]*dρ[2]*d7x[1] + 28*dx[1]^6*d2x[1]*dρ[7] + dρ[1]*d8x[1] +
dx[1]^8*dρ[8])

yH(x::AbstractVector{Float64}, ξ::Float64, ρ::Float64) = sin(x[2])*sin(ξ)*ρ

dyH1(x::AbstractVector{Float64}, dx::AbstractVector{Float64}, ξ::Float64, dξ::AbstractVector{Float64}, ρ::Float64, dρ::AbstractVector{Float64}) = cos(x[2])*sin(ξ)*ρ*dx[2] + cos(ξ)*sin(x[2])*ρ*dξ[1] +
sin(x[2])*sin(ξ)*dx[1]*dρ[1]

dyH2(x::AbstractVector{Float64}, dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, ξ::Float64, dξ::AbstractVector{Float64}, ρ::Float64, dρ::AbstractVector{Float64}) = 2*dx[1]*(cos(x[2])*sin(ξ)*dx[2] +
cos(ξ)*sin(x[2])*dξ[1])*dρ[1] + ρ*(2*cos(x[2])*cos(ξ)*dx[2]*dξ[1] + sin(ξ)*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) + sin(x[2])*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2])) +
sin(x[2])*sin(ξ)*(dρ[1]*d2x[1] + dx[1]^2*dρ[2])

dyH3(x::AbstractVector{Float64}, dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, ξ::Float64, dξ::AbstractVector{Float64}, ρ::Float64, dρ::AbstractVector{Float64}) = 3*dx[1]*dρ[1]*(2*cos(x[2])*cos(ξ)*dx[2]*dξ[1] + sin(ξ)*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) + sin(x[2])*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2])) + 3*(cos(x[2])*sin(ξ)*dx[2] +
cos(ξ)*sin(x[2])*dξ[1])*(dρ[1]*d2x[1] + dx[1]^2*dρ[2]) + ρ*(3*cos(ξ)*dξ[1]*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) + 3*cos(x[2])*dx[2]*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2]) +
sin(ξ)*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) + sin(x[2])*(-(cos(ξ)*dξ[1]^3) - 3*sin(ξ)*dξ[1]*dξ[2] + cos(ξ)*dξ[3])) +
sin(x[2])*sin(ξ)*(3*dx[1]*d2x[1]*dρ[2] + dρ[1]*d3x[1] + dx[1]^3*dρ[3])

dyH4(x::AbstractVector{Float64}, dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, d4x::AbstractVector{Float64}, ξ::Float64, dξ::AbstractVector{Float64}, ρ::Float64, dρ::AbstractVector{Float64}) = 6*(2*cos(x[2])*cos(ξ)*dx[2]*dξ[1] + sin(ξ)*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) + sin(x[2])*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2]))*(dρ[1]*d2x[1] + dx[1]^2*dρ[2]) +
4*dx[1]*dρ[1]*(3*cos(ξ)*dξ[1]*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) + 3*cos(x[2])*dx[2]*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2]) + sin(ξ)*(-(cos(x[2])*dx[2]^3) -
3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) + sin(x[2])*(-(cos(ξ)*dξ[1]^3) - 3*sin(ξ)*dξ[1]*dξ[2] + cos(ξ)*dξ[3])) + 4*(cos(x[2])*sin(ξ)*dx[2] +
cos(ξ)*sin(x[2])*dξ[1])*(3*dx[1]*d2x[1]*dρ[2] + dρ[1]*d3x[1] + dx[1]^3*dρ[3]) + ρ*(6*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2]) +
4*cos(ξ)*dξ[1]*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) + 4*cos(x[2])*dx[2]*(-(cos(ξ)*dξ[1]^3) - 3*sin(ξ)*dξ[1]*dξ[2] + cos(ξ)*dξ[3]) +
sin(ξ)*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + sin(x[2])*(sin(ξ)*dξ[1]^4 -
6*cos(ξ)*dξ[1]^2*dξ[2] - 3*sin(ξ)*dξ[2]^2 - 4*sin(ξ)*dξ[1]*dξ[3] + cos(ξ)*dξ[4])) + sin(x[2])*sin(ξ)*(3*d2x[1]^2*dρ[2] + 4*dx[1]*dρ[2]*d3x[1] + 6*dx[1]^2*d2x[1]*dρ[3] +
dρ[1]*d4x[1] + dx[1]^4*dρ[4])

dyH5(x::AbstractVector{Float64}, dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, d4x::AbstractVector{Float64}, d5x::AbstractVector{Float64}, ξ::Float64, dξ::AbstractVector{Float64}, ρ::Float64, dρ::AbstractVector{Float64}) = 10*(dρ[1]*d2x[1] + dx[1]^2*dρ[2])*(3*cos(ξ)*dξ[1]*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) + 3*cos(x[2])*dx[2]*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2]) +
sin(ξ)*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) + sin(x[2])*(-(cos(ξ)*dξ[1]^3) - 3*sin(ξ)*dξ[1]*dξ[2] + cos(ξ)*dξ[3])) +
10*(2*cos(x[2])*cos(ξ)*dx[2]*dξ[1] + sin(ξ)*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) + sin(x[2])*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2]))*(3*dx[1]*d2x[1]*dρ[2] + dρ[1]*d3x[1] +
dx[1]^3*dρ[3]) + 5*dx[1]*dρ[1]*(6*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2]) + 4*cos(ξ)*dξ[1]*(-(cos(x[2])*dx[2]^3) -
3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) + 4*cos(x[2])*dx[2]*(-(cos(ξ)*dξ[1]^3) - 3*sin(ξ)*dξ[1]*dξ[2] + cos(ξ)*dξ[3]) + sin(ξ)*(sin(x[2])*dx[2]^4 -
6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + sin(x[2])*(sin(ξ)*dξ[1]^4 - 6*cos(ξ)*dξ[1]^2*dξ[2] - 3*sin(ξ)*dξ[2]^2 -
4*sin(ξ)*dξ[1]*dξ[3] + cos(ξ)*dξ[4])) + 5*(cos(x[2])*sin(ξ)*dx[2] + cos(ξ)*sin(x[2])*dξ[1])*(3*d2x[1]^2*dρ[2] + 4*dx[1]*dρ[2]*d3x[1] + 6*dx[1]^2*d2x[1]*dρ[3] + dρ[1]*d4x[1] +
dx[1]^4*dρ[4]) + ρ*(10*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2])*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) + 10*(-(sin(x[2])*dx[2]^2) +
cos(x[2])*d2x[2])*(-(cos(ξ)*dξ[1]^3) - 3*sin(ξ)*dξ[1]*dξ[2] + cos(ξ)*dξ[3]) + 5*cos(ξ)*dξ[1]*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 -
4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + 5*cos(x[2])*dx[2]*(sin(ξ)*dξ[1]^4 - 6*cos(ξ)*dξ[1]^2*dξ[2] - 3*sin(ξ)*dξ[2]^2 - 4*sin(ξ)*dξ[1]*dξ[3] + cos(ξ)*dξ[4]) +
sin(ξ)*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] - 15*cos(x[2])*dx[2]*d2x[2]^2 - 10*cos(x[2])*dx[2]^2*d3x[2] - 10*sin(x[2])*d2x[2]*d3x[2] - 5*sin(x[2])*dx[2]*d4x[2] +
cos(x[2])*d5x[2]) + sin(x[2])*(cos(ξ)*dξ[1]^5 + 10*sin(ξ)*dξ[1]^3*dξ[2] - 15*cos(ξ)*dξ[1]*dξ[2]^2 - 10*cos(ξ)*dξ[1]^2*dξ[3] - 10*sin(ξ)*dξ[2]*dξ[3] - 5*sin(ξ)*dξ[1]*dξ[4] +
cos(ξ)*dξ[5])) + sin(x[2])*sin(ξ)*(10*d2x[1]*dρ[2]*d3x[1] + 15*dx[1]*d2x[1]^2*dρ[3] + 10*dx[1]^2*d3x[1]*dρ[3] + 5*dx[1]*dρ[2]*d4x[1] + 10*dx[1]^3*d2x[1]*dρ[4] + dρ[1]*d5x[1] +
dx[1]^5*dρ[5])

dyH6(x::AbstractVector{Float64}, dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, d4x::AbstractVector{Float64}, d5x::AbstractVector{Float64}, d6x::AbstractVector{Float64}, ξ::Float64, dξ::AbstractVector{Float64}, ρ::Float64, dρ::AbstractVector{Float64}) = 20*(3*cos(ξ)*dξ[1]*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) + 3*cos(x[2])*dx[2]*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2]) +
sin(ξ)*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) + sin(x[2])*(-(cos(ξ)*dξ[1]^3) - 3*sin(ξ)*dξ[1]*dξ[2] + cos(ξ)*dξ[3]))*(3*dx[1]*d2x[1]*dρ[2] +
dρ[1]*d3x[1] + dx[1]^3*dρ[3]) + 15*(dρ[1]*d2x[1] + dx[1]^2*dρ[2])*(6*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2]) +
4*cos(ξ)*dξ[1]*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) + 4*cos(x[2])*dx[2]*(-(cos(ξ)*dξ[1]^3) - 3*sin(ξ)*dξ[1]*dξ[2] + cos(ξ)*dξ[3]) +
sin(ξ)*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + sin(x[2])*(sin(ξ)*dξ[1]^4 -
6*cos(ξ)*dξ[1]^2*dξ[2] - 3*sin(ξ)*dξ[2]^2 - 4*sin(ξ)*dξ[1]*dξ[3] + cos(ξ)*dξ[4])) + 15*(2*cos(x[2])*cos(ξ)*dx[2]*dξ[1] + sin(ξ)*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) +
sin(x[2])*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2]))*(3*d2x[1]^2*dρ[2] + 4*dx[1]*dρ[2]*d3x[1] + 6*dx[1]^2*d2x[1]*dρ[3] + dρ[1]*d4x[1] + dx[1]^4*dρ[4]) +
6*dx[1]*dρ[1]*(10*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2])*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) + 10*(-(sin(x[2])*dx[2]^2) +
cos(x[2])*d2x[2])*(-(cos(ξ)*dξ[1]^3) - 3*sin(ξ)*dξ[1]*dξ[2] + cos(ξ)*dξ[3]) + 5*cos(ξ)*dξ[1]*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 -
4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + 5*cos(x[2])*dx[2]*(sin(ξ)*dξ[1]^4 - 6*cos(ξ)*dξ[1]^2*dξ[2] - 3*sin(ξ)*dξ[2]^2 - 4*sin(ξ)*dξ[1]*dξ[3] + cos(ξ)*dξ[4]) +
sin(ξ)*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] - 15*cos(x[2])*dx[2]*d2x[2]^2 - 10*cos(x[2])*dx[2]^2*d3x[2] - 10*sin(x[2])*d2x[2]*d3x[2] - 5*sin(x[2])*dx[2]*d4x[2] +
cos(x[2])*d5x[2]) + sin(x[2])*(cos(ξ)*dξ[1]^5 + 10*sin(ξ)*dξ[1]^3*dξ[2] - 15*cos(ξ)*dξ[1]*dξ[2]^2 - 10*cos(ξ)*dξ[1]^2*dξ[3] - 10*sin(ξ)*dξ[2]*dξ[3] - 5*sin(ξ)*dξ[1]*dξ[4] +
cos(ξ)*dξ[5])) + 6*(cos(x[2])*sin(ξ)*dx[2] + cos(ξ)*sin(x[2])*dξ[1])*(10*d2x[1]*dρ[2]*d3x[1] + 15*dx[1]*d2x[1]^2*dρ[3] + 10*dx[1]^2*d3x[1]*dρ[3] + 5*dx[1]*dρ[2]*d4x[1] +
10*dx[1]^3*d2x[1]*dρ[4] + dρ[1]*d5x[1] + dx[1]^5*dρ[5]) + ρ*(20*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2])*(-(cos(ξ)*dξ[1]^3) - 3*sin(ξ)*dξ[1]*dξ[2] +
cos(ξ)*dξ[3]) + 15*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2])*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) +
15*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(sin(ξ)*dξ[1]^4 - 6*cos(ξ)*dξ[1]^2*dξ[2] - 3*sin(ξ)*dξ[2]^2 - 4*sin(ξ)*dξ[1]*dξ[3] + cos(ξ)*dξ[4]) +
6*cos(ξ)*dξ[1]*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] - 15*cos(x[2])*dx[2]*d2x[2]^2 - 10*cos(x[2])*dx[2]^2*d3x[2] - 10*sin(x[2])*d2x[2]*d3x[2] -
5*sin(x[2])*dx[2]*d4x[2] + cos(x[2])*d5x[2]) + 6*cos(x[2])*dx[2]*(cos(ξ)*dξ[1]^5 + 10*sin(ξ)*dξ[1]^3*dξ[2] - 15*cos(ξ)*dξ[1]*dξ[2]^2 - 10*cos(ξ)*dξ[1]^2*dξ[3] -
10*sin(ξ)*dξ[2]*dξ[3] - 5*sin(ξ)*dξ[1]*dξ[4] + cos(ξ)*dξ[5]) + sin(ξ)*(-(sin(x[2])*dx[2]^6) + 15*cos(x[2])*dx[2]^4*d2x[2] + 45*sin(x[2])*dx[2]^2*d2x[2]^2 -
15*cos(x[2])*d2x[2]^3 + 20*sin(x[2])*dx[2]^3*d3x[2] - 60*cos(x[2])*dx[2]*d2x[2]*d3x[2] - 10*sin(x[2])*d3x[2]^2 - 15*cos(x[2])*dx[2]^2*d4x[2] - 15*sin(x[2])*d2x[2]*d4x[2] -
6*sin(x[2])*dx[2]*d5x[2] + cos(x[2])*d6x[2]) + sin(x[2])*(-(sin(ξ)*dξ[1]^6) + 15*cos(ξ)*dξ[1]^4*dξ[2] + 45*sin(ξ)*dξ[1]^2*dξ[2]^2 - 15*cos(ξ)*dξ[2]^3 +
20*sin(ξ)*dξ[1]^3*dξ[3] - 60*cos(ξ)*dξ[1]*dξ[2]*dξ[3] - 10*sin(ξ)*dξ[3]^2 - 15*cos(ξ)*dξ[1]^2*dξ[4] - 15*sin(ξ)*dξ[2]*dξ[4] - 6*sin(ξ)*dξ[1]*dξ[5] + cos(ξ)*dξ[6])) +
sin(x[2])*sin(ξ)*(10*dρ[2]*d3x[1]^2 + 15*d2x[1]^3*dρ[3] + 60*dx[1]*d2x[1]*d3x[1]*dρ[3] + 15*d2x[1]*dρ[2]*d4x[1] + 15*dx[1]^2*dρ[3]*d4x[1] + 45*dx[1]^2*d2x[1]^2*dρ[4] +
20*dx[1]^3*d3x[1]*dρ[4] + 6*dx[1]*dρ[2]*d5x[1] + 15*dx[1]^4*d2x[1]*dρ[5] + dρ[1]*d6x[1] + dx[1]^6*dρ[6])

dyH7(x::AbstractVector{Float64}, dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, d4x::AbstractVector{Float64}, d5x::AbstractVector{Float64}, d6x::AbstractVector{Float64}, d7x::AbstractVector{Float64}, ξ::Float64, dξ::AbstractVector{Float64}, ρ::Float64, dρ::AbstractVector{Float64}) = 35*(3*dx[1]*d2x[1]*dρ[2] + dρ[1]*d3x[1] + dx[1]^3*dρ[3])*(6*(-(sin(x[2])*dx[2]^2) +
cos(x[2])*d2x[2])*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2]) + 4*cos(ξ)*dξ[1]*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) +
4*cos(x[2])*dx[2]*(-(cos(ξ)*dξ[1]^3) - 3*sin(ξ)*dξ[1]*dξ[2] + cos(ξ)*dξ[3]) + sin(ξ)*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 -
4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + sin(x[2])*(sin(ξ)*dξ[1]^4 - 6*cos(ξ)*dξ[1]^2*dξ[2] - 3*sin(ξ)*dξ[2]^2 - 4*sin(ξ)*dξ[1]*dξ[3] + cos(ξ)*dξ[4])) +
35*(3*cos(ξ)*dξ[1]*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) + 3*cos(x[2])*dx[2]*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2]) + sin(ξ)*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] +
cos(x[2])*d3x[2]) + sin(x[2])*(-(cos(ξ)*dξ[1]^3) - 3*sin(ξ)*dξ[1]*dξ[2] + cos(ξ)*dξ[3]))*(3*d2x[1]^2*dρ[2] + 4*dx[1]*dρ[2]*d3x[1] + 6*dx[1]^2*d2x[1]*dρ[3] + dρ[1]*d4x[1] +
dx[1]^4*dρ[4]) + 21*(dρ[1]*d2x[1] + dx[1]^2*dρ[2])*(10*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2])*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) +
10*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(-(cos(ξ)*dξ[1]^3) - 3*sin(ξ)*dξ[1]*dξ[2] + cos(ξ)*dξ[3]) + 5*cos(ξ)*dξ[1]*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] -
3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + 5*cos(x[2])*dx[2]*(sin(ξ)*dξ[1]^4 - 6*cos(ξ)*dξ[1]^2*dξ[2] - 3*sin(ξ)*dξ[2]^2 - 4*sin(ξ)*dξ[1]*dξ[3] +
cos(ξ)*dξ[4]) + sin(ξ)*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] - 15*cos(x[2])*dx[2]*d2x[2]^2 - 10*cos(x[2])*dx[2]^2*d3x[2] - 10*sin(x[2])*d2x[2]*d3x[2] -
5*sin(x[2])*dx[2]*d4x[2] + cos(x[2])*d5x[2]) + sin(x[2])*(cos(ξ)*dξ[1]^5 + 10*sin(ξ)*dξ[1]^3*dξ[2] - 15*cos(ξ)*dξ[1]*dξ[2]^2 - 10*cos(ξ)*dξ[1]^2*dξ[3] - 10*sin(ξ)*dξ[2]*dξ[3] -
5*sin(ξ)*dξ[1]*dξ[4] + cos(ξ)*dξ[5])) + 21*(2*cos(x[2])*cos(ξ)*dx[2]*dξ[1] + sin(ξ)*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) +
sin(x[2])*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2]))*(10*d2x[1]*dρ[2]*d3x[1] + 15*dx[1]*d2x[1]^2*dρ[3] + 10*dx[1]^2*d3x[1]*dρ[3] + 5*dx[1]*dρ[2]*d4x[1] + 10*dx[1]^3*d2x[1]*dρ[4] +
dρ[1]*d5x[1] + dx[1]^5*dρ[5]) + 7*dx[1]*dρ[1]*(20*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2])*(-(cos(ξ)*dξ[1]^3) - 3*sin(ξ)*dξ[1]*dξ[2] +
cos(ξ)*dξ[3]) + 15*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2])*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) +
15*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(sin(ξ)*dξ[1]^4 - 6*cos(ξ)*dξ[1]^2*dξ[2] - 3*sin(ξ)*dξ[2]^2 - 4*sin(ξ)*dξ[1]*dξ[3] + cos(ξ)*dξ[4]) +
6*cos(ξ)*dξ[1]*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] - 15*cos(x[2])*dx[2]*d2x[2]^2 - 10*cos(x[2])*dx[2]^2*d3x[2] - 10*sin(x[2])*d2x[2]*d3x[2] -
5*sin(x[2])*dx[2]*d4x[2] + cos(x[2])*d5x[2]) + 6*cos(x[2])*dx[2]*(cos(ξ)*dξ[1]^5 + 10*sin(ξ)*dξ[1]^3*dξ[2] - 15*cos(ξ)*dξ[1]*dξ[2]^2 - 10*cos(ξ)*dξ[1]^2*dξ[3] -
10*sin(ξ)*dξ[2]*dξ[3] - 5*sin(ξ)*dξ[1]*dξ[4] + cos(ξ)*dξ[5]) + sin(ξ)*(-(sin(x[2])*dx[2]^6) + 15*cos(x[2])*dx[2]^4*d2x[2] + 45*sin(x[2])*dx[2]^2*d2x[2]^2 -
15*cos(x[2])*d2x[2]^3 + 20*sin(x[2])*dx[2]^3*d3x[2] - 60*cos(x[2])*dx[2]*d2x[2]*d3x[2] - 10*sin(x[2])*d3x[2]^2 - 15*cos(x[2])*dx[2]^2*d4x[2] - 15*sin(x[2])*d2x[2]*d4x[2] -
6*sin(x[2])*dx[2]*d5x[2] + cos(x[2])*d6x[2]) + sin(x[2])*(-(sin(ξ)*dξ[1]^6) + 15*cos(ξ)*dξ[1]^4*dξ[2] + 45*sin(ξ)*dξ[1]^2*dξ[2]^2 - 15*cos(ξ)*dξ[2]^3 +
20*sin(ξ)*dξ[1]^3*dξ[3] - 60*cos(ξ)*dξ[1]*dξ[2]*dξ[3] - 10*sin(ξ)*dξ[3]^2 - 15*cos(ξ)*dξ[1]^2*dξ[4] - 15*sin(ξ)*dξ[2]*dξ[4] - 6*sin(ξ)*dξ[1]*dξ[5] + cos(ξ)*dξ[6])) +
7*(cos(x[2])*sin(ξ)*dx[2] + cos(ξ)*sin(x[2])*dξ[1])*(10*dρ[2]*d3x[1]^2 + 15*d2x[1]^3*dρ[3] + 60*dx[1]*d2x[1]*d3x[1]*dρ[3] + 15*d2x[1]*dρ[2]*d4x[1] + 15*dx[1]^2*dρ[3]*d4x[1] +
45*dx[1]^2*d2x[1]^2*dρ[4] + 20*dx[1]^3*d3x[1]*dρ[4] + 6*dx[1]*dρ[2]*d5x[1] + 15*dx[1]^4*d2x[1]*dρ[5] + dρ[1]*d6x[1] + dx[1]^6*dρ[6]) + ρ*(35*(-(cos(ξ)*dξ[1]^3) -
3*sin(ξ)*dξ[1]*dξ[2] + cos(ξ)*dξ[3])*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) +
35*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2])*(sin(ξ)*dξ[1]^4 - 6*cos(ξ)*dξ[1]^2*dξ[2] - 3*sin(ξ)*dξ[2]^2 - 4*sin(ξ)*dξ[1]*dξ[3] + cos(ξ)*dξ[4]) +
21*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2])*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] - 15*cos(x[2])*dx[2]*d2x[2]^2 - 10*cos(x[2])*dx[2]^2*d3x[2] -
10*sin(x[2])*d2x[2]*d3x[2] - 5*sin(x[2])*dx[2]*d4x[2] + cos(x[2])*d5x[2]) + 21*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(cos(ξ)*dξ[1]^5 + 10*sin(ξ)*dξ[1]^3*dξ[2] -
15*cos(ξ)*dξ[1]*dξ[2]^2 - 10*cos(ξ)*dξ[1]^2*dξ[3] - 10*sin(ξ)*dξ[2]*dξ[3] - 5*sin(ξ)*dξ[1]*dξ[4] + cos(ξ)*dξ[5]) + 7*cos(ξ)*dξ[1]*(-(sin(x[2])*dx[2]^6) +
15*cos(x[2])*dx[2]^4*d2x[2] + 45*sin(x[2])*dx[2]^2*d2x[2]^2 - 15*cos(x[2])*d2x[2]^3 + 20*sin(x[2])*dx[2]^3*d3x[2] - 60*cos(x[2])*dx[2]*d2x[2]*d3x[2] - 10*sin(x[2])*d3x[2]^2 -
15*cos(x[2])*dx[2]^2*d4x[2] - 15*sin(x[2])*d2x[2]*d4x[2] - 6*sin(x[2])*dx[2]*d5x[2] + cos(x[2])*d6x[2]) + 7*cos(x[2])*dx[2]*(-(sin(ξ)*dξ[1]^6) + 15*cos(ξ)*dξ[1]^4*dξ[2] +
45*sin(ξ)*dξ[1]^2*dξ[2]^2 - 15*cos(ξ)*dξ[2]^3 + 20*sin(ξ)*dξ[1]^3*dξ[3] - 60*cos(ξ)*dξ[1]*dξ[2]*dξ[3] - 10*sin(ξ)*dξ[3]^2 - 15*cos(ξ)*dξ[1]^2*dξ[4] - 15*sin(ξ)*dξ[2]*dξ[4] -
6*sin(ξ)*dξ[1]*dξ[5] + cos(ξ)*dξ[6]) + sin(ξ)*(-(cos(x[2])*dx[2]^7) - 21*sin(x[2])*dx[2]^5*d2x[2] + 105*cos(x[2])*dx[2]^3*d2x[2]^2 + 105*sin(x[2])*dx[2]*d2x[2]^3 +
35*cos(x[2])*dx[2]^4*d3x[2] + 210*sin(x[2])*dx[2]^2*d2x[2]*d3x[2] - 105*cos(x[2])*d2x[2]^2*d3x[2] - 70*cos(x[2])*dx[2]*d3x[2]^2 + 35*sin(x[2])*dx[2]^3*d4x[2] -
105*cos(x[2])*dx[2]*d2x[2]*d4x[2] - 35*sin(x[2])*d3x[2]*d4x[2] - 21*cos(x[2])*dx[2]^2*d5x[2] - 21*sin(x[2])*d2x[2]*d5x[2] - 7*sin(x[2])*dx[2]*d6x[2] + cos(x[2])*d7x[2]) +
sin(x[2])*(-(cos(ξ)*dξ[1]^7) - 21*sin(ξ)*dξ[1]^5*dξ[2] + 105*cos(ξ)*dξ[1]^3*dξ[2]^2 + 105*sin(ξ)*dξ[1]*dξ[2]^3 + 35*cos(ξ)*dξ[1]^4*dξ[3] + 210*sin(ξ)*dξ[1]^2*dξ[2]*dξ[3] -
105*cos(ξ)*dξ[2]^2*dξ[3] - 70*cos(ξ)*dξ[1]*dξ[3]^2 + 35*sin(ξ)*dξ[1]^3*dξ[4] - 105*cos(ξ)*dξ[1]*dξ[2]*dξ[4] - 35*sin(ξ)*dξ[3]*dξ[4] - 21*cos(ξ)*dξ[1]^2*dξ[5] -
21*sin(ξ)*dξ[2]*dξ[5] - 7*sin(ξ)*dξ[1]*dξ[6] + cos(ξ)*dξ[7])) + sin(x[2])*sin(ξ)*(105*d2x[1]^2*d3x[1]*dρ[3] + 70*dx[1]*d3x[1]^2*dρ[3] + 35*dρ[2]*d3x[1]*d4x[1] +
105*dx[1]*d2x[1]*dρ[3]*d4x[1] + 105*dx[1]*d2x[1]^3*dρ[4] + 210*dx[1]^2*d2x[1]*d3x[1]*dρ[4] + 35*dx[1]^3*d4x[1]*dρ[4] + 21*d2x[1]*dρ[2]*d5x[1] + 21*dx[1]^2*dρ[3]*d5x[1] +
105*dx[1]^3*d2x[1]^2*dρ[5] + 35*dx[1]^4*d3x[1]*dρ[5] + 7*dx[1]*dρ[2]*d6x[1] + 21*dx[1]^5*d2x[1]*dρ[6] + dρ[1]*d7x[1] + dx[1]^7*dρ[7])

dyH8(x::AbstractVector{Float64}, dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, d4x::AbstractVector{Float64}, d5x::AbstractVector{Float64}, d6x::AbstractVector{Float64}, d7x::AbstractVector{Float64}, d8x::AbstractVector{Float64}, ξ::Float64, dξ::AbstractVector{Float64}, ρ::Float64, dρ::AbstractVector{Float64}) = 70*(6*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2]) +
4*cos(ξ)*dξ[1]*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) + 4*cos(x[2])*dx[2]*(-(cos(ξ)*dξ[1]^3) - 3*sin(ξ)*dξ[1]*dξ[2] + cos(ξ)*dξ[3]) +
sin(ξ)*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + sin(x[2])*(sin(ξ)*dξ[1]^4 -
6*cos(ξ)*dξ[1]^2*dξ[2] - 3*sin(ξ)*dξ[2]^2 - 4*sin(ξ)*dξ[1]*dξ[3] + cos(ξ)*dξ[4]))*(3*d2x[1]^2*dρ[2] + 4*dx[1]*dρ[2]*d3x[1] + 6*dx[1]^2*d2x[1]*dρ[3] + dρ[1]*d4x[1] +
dx[1]^4*dρ[4]) + 56*(3*dx[1]*d2x[1]*dρ[2] + dρ[1]*d3x[1] + dx[1]^3*dρ[3])*(10*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2])*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] +
cos(x[2])*d3x[2]) + 10*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(-(cos(ξ)*dξ[1]^3) - 3*sin(ξ)*dξ[1]*dξ[2] + cos(ξ)*dξ[3]) + 5*cos(ξ)*dξ[1]*(sin(x[2])*dx[2]^4 -
6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + 5*cos(x[2])*dx[2]*(sin(ξ)*dξ[1]^4 - 6*cos(ξ)*dξ[1]^2*dξ[2] -
3*sin(ξ)*dξ[2]^2 - 4*sin(ξ)*dξ[1]*dξ[3] + cos(ξ)*dξ[4]) + sin(ξ)*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] - 15*cos(x[2])*dx[2]*d2x[2]^2 - 10*cos(x[2])*dx[2]^2*d3x[2] -
10*sin(x[2])*d2x[2]*d3x[2] - 5*sin(x[2])*dx[2]*d4x[2] + cos(x[2])*d5x[2]) + sin(x[2])*(cos(ξ)*dξ[1]^5 + 10*sin(ξ)*dξ[1]^3*dξ[2] - 15*cos(ξ)*dξ[1]*dξ[2]^2 -
10*cos(ξ)*dξ[1]^2*dξ[3] - 10*sin(ξ)*dξ[2]*dξ[3] - 5*sin(ξ)*dξ[1]*dξ[4] + cos(ξ)*dξ[5])) + 56*(3*cos(ξ)*dξ[1]*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) +
3*cos(x[2])*dx[2]*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2]) + sin(ξ)*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) + sin(x[2])*(-(cos(ξ)*dξ[1]^3) -
3*sin(ξ)*dξ[1]*dξ[2] + cos(ξ)*dξ[3]))*(10*d2x[1]*dρ[2]*d3x[1] + 15*dx[1]*d2x[1]^2*dρ[3] + 10*dx[1]^2*d3x[1]*dρ[3] + 5*dx[1]*dρ[2]*d4x[1] + 10*dx[1]^3*d2x[1]*dρ[4] +
dρ[1]*d5x[1] + dx[1]^5*dρ[5]) + 28*(dρ[1]*d2x[1] + dx[1]^2*dρ[2])*(20*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2])*(-(cos(ξ)*dξ[1]^3) -
3*sin(ξ)*dξ[1]*dξ[2] + cos(ξ)*dξ[3]) + 15*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2])*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] +
cos(x[2])*d4x[2]) + 15*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(sin(ξ)*dξ[1]^4 - 6*cos(ξ)*dξ[1]^2*dξ[2] - 3*sin(ξ)*dξ[2]^2 - 4*sin(ξ)*dξ[1]*dξ[3] + cos(ξ)*dξ[4]) +
6*cos(ξ)*dξ[1]*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] - 15*cos(x[2])*dx[2]*d2x[2]^2 - 10*cos(x[2])*dx[2]^2*d3x[2] - 10*sin(x[2])*d2x[2]*d3x[2] -
5*sin(x[2])*dx[2]*d4x[2] + cos(x[2])*d5x[2]) + 6*cos(x[2])*dx[2]*(cos(ξ)*dξ[1]^5 + 10*sin(ξ)*dξ[1]^3*dξ[2] - 15*cos(ξ)*dξ[1]*dξ[2]^2 - 10*cos(ξ)*dξ[1]^2*dξ[3] -
10*sin(ξ)*dξ[2]*dξ[3] - 5*sin(ξ)*dξ[1]*dξ[4] + cos(ξ)*dξ[5]) + sin(ξ)*(-(sin(x[2])*dx[2]^6) + 15*cos(x[2])*dx[2]^4*d2x[2] + 45*sin(x[2])*dx[2]^2*d2x[2]^2 -
15*cos(x[2])*d2x[2]^3 + 20*sin(x[2])*dx[2]^3*d3x[2] - 60*cos(x[2])*dx[2]*d2x[2]*d3x[2] - 10*sin(x[2])*d3x[2]^2 - 15*cos(x[2])*dx[2]^2*d4x[2] - 15*sin(x[2])*d2x[2]*d4x[2] -
6*sin(x[2])*dx[2]*d5x[2] + cos(x[2])*d6x[2]) + sin(x[2])*(-(sin(ξ)*dξ[1]^6) + 15*cos(ξ)*dξ[1]^4*dξ[2] + 45*sin(ξ)*dξ[1]^2*dξ[2]^2 - 15*cos(ξ)*dξ[2]^3 +
20*sin(ξ)*dξ[1]^3*dξ[3] - 60*cos(ξ)*dξ[1]*dξ[2]*dξ[3] - 10*sin(ξ)*dξ[3]^2 - 15*cos(ξ)*dξ[1]^2*dξ[4] - 15*sin(ξ)*dξ[2]*dξ[4] - 6*sin(ξ)*dξ[1]*dξ[5] + cos(ξ)*dξ[6])) +
28*(2*cos(x[2])*cos(ξ)*dx[2]*dξ[1] + sin(ξ)*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) + sin(x[2])*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2]))*(10*dρ[2]*d3x[1]^2 + 15*d2x[1]^3*dρ[3] +
60*dx[1]*d2x[1]*d3x[1]*dρ[3] + 15*d2x[1]*dρ[2]*d4x[1] + 15*dx[1]^2*dρ[3]*d4x[1] + 45*dx[1]^2*d2x[1]^2*dρ[4] + 20*dx[1]^3*d3x[1]*dρ[4] + 6*dx[1]*dρ[2]*d5x[1] +
15*dx[1]^4*d2x[1]*dρ[5] + dρ[1]*d6x[1] + dx[1]^6*dρ[6]) + 8*dx[1]*dρ[1]*(35*(-(cos(ξ)*dξ[1]^3) - 3*sin(ξ)*dξ[1]*dξ[2] + cos(ξ)*dξ[3])*(sin(x[2])*dx[2]^4 -
6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + 35*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] +
cos(x[2])*d3x[2])*(sin(ξ)*dξ[1]^4 - 6*cos(ξ)*dξ[1]^2*dξ[2] - 3*sin(ξ)*dξ[2]^2 - 4*sin(ξ)*dξ[1]*dξ[3] + cos(ξ)*dξ[4]) + 21*(-(sin(ξ)*dξ[1]^2) +
cos(ξ)*dξ[2])*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] - 15*cos(x[2])*dx[2]*d2x[2]^2 - 10*cos(x[2])*dx[2]^2*d3x[2] - 10*sin(x[2])*d2x[2]*d3x[2] -
5*sin(x[2])*dx[2]*d4x[2] + cos(x[2])*d5x[2]) + 21*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(cos(ξ)*dξ[1]^5 + 10*sin(ξ)*dξ[1]^3*dξ[2] - 15*cos(ξ)*dξ[1]*dξ[2]^2 -
10*cos(ξ)*dξ[1]^2*dξ[3] - 10*sin(ξ)*dξ[2]*dξ[3] - 5*sin(ξ)*dξ[1]*dξ[4] + cos(ξ)*dξ[5]) + 7*cos(ξ)*dξ[1]*(-(sin(x[2])*dx[2]^6) + 15*cos(x[2])*dx[2]^4*d2x[2] +
45*sin(x[2])*dx[2]^2*d2x[2]^2 - 15*cos(x[2])*d2x[2]^3 + 20*sin(x[2])*dx[2]^3*d3x[2] - 60*cos(x[2])*dx[2]*d2x[2]*d3x[2] - 10*sin(x[2])*d3x[2]^2 - 15*cos(x[2])*dx[2]^2*d4x[2] -
15*sin(x[2])*d2x[2]*d4x[2] - 6*sin(x[2])*dx[2]*d5x[2] + cos(x[2])*d6x[2]) + 7*cos(x[2])*dx[2]*(-(sin(ξ)*dξ[1]^6) + 15*cos(ξ)*dξ[1]^4*dξ[2] + 45*sin(ξ)*dξ[1]^2*dξ[2]^2 -
15*cos(ξ)*dξ[2]^3 + 20*sin(ξ)*dξ[1]^3*dξ[3] - 60*cos(ξ)*dξ[1]*dξ[2]*dξ[3] - 10*sin(ξ)*dξ[3]^2 - 15*cos(ξ)*dξ[1]^2*dξ[4] - 15*sin(ξ)*dξ[2]*dξ[4] - 6*sin(ξ)*dξ[1]*dξ[5] +
cos(ξ)*dξ[6]) + sin(ξ)*(-(cos(x[2])*dx[2]^7) - 21*sin(x[2])*dx[2]^5*d2x[2] + 105*cos(x[2])*dx[2]^3*d2x[2]^2 + 105*sin(x[2])*dx[2]*d2x[2]^3 + 35*cos(x[2])*dx[2]^4*d3x[2] +
210*sin(x[2])*dx[2]^2*d2x[2]*d3x[2] - 105*cos(x[2])*d2x[2]^2*d3x[2] - 70*cos(x[2])*dx[2]*d3x[2]^2 + 35*sin(x[2])*dx[2]^3*d4x[2] - 105*cos(x[2])*dx[2]*d2x[2]*d4x[2] -
35*sin(x[2])*d3x[2]*d4x[2] - 21*cos(x[2])*dx[2]^2*d5x[2] - 21*sin(x[2])*d2x[2]*d5x[2] - 7*sin(x[2])*dx[2]*d6x[2] + cos(x[2])*d7x[2]) + sin(x[2])*(-(cos(ξ)*dξ[1]^7) -
21*sin(ξ)*dξ[1]^5*dξ[2] + 105*cos(ξ)*dξ[1]^3*dξ[2]^2 + 105*sin(ξ)*dξ[1]*dξ[2]^3 + 35*cos(ξ)*dξ[1]^4*dξ[3] + 210*sin(ξ)*dξ[1]^2*dξ[2]*dξ[3] - 105*cos(ξ)*dξ[2]^2*dξ[3] -
70*cos(ξ)*dξ[1]*dξ[3]^2 + 35*sin(ξ)*dξ[1]^3*dξ[4] - 105*cos(ξ)*dξ[1]*dξ[2]*dξ[4] - 35*sin(ξ)*dξ[3]*dξ[4] - 21*cos(ξ)*dξ[1]^2*dξ[5] - 21*sin(ξ)*dξ[2]*dξ[5] -
7*sin(ξ)*dξ[1]*dξ[6] + cos(ξ)*dξ[7])) + 8*(cos(x[2])*sin(ξ)*dx[2] + cos(ξ)*sin(x[2])*dξ[1])*(105*d2x[1]^2*d3x[1]*dρ[3] + 70*dx[1]*d3x[1]^2*dρ[3] + 35*dρ[2]*d3x[1]*d4x[1] +
105*dx[1]*d2x[1]*dρ[3]*d4x[1] + 105*dx[1]*d2x[1]^3*dρ[4] + 210*dx[1]^2*d2x[1]*d3x[1]*dρ[4] + 35*dx[1]^3*d4x[1]*dρ[4] + 21*d2x[1]*dρ[2]*d5x[1] + 21*dx[1]^2*dρ[3]*d5x[1] +
105*dx[1]^3*d2x[1]^2*dρ[5] + 35*dx[1]^4*d3x[1]*dρ[5] + 7*dx[1]*dρ[2]*d6x[1] + 21*dx[1]^5*d2x[1]*dρ[6] + dρ[1]*d7x[1] + dx[1]^7*dρ[7]) + ρ*(70*(sin(x[2])*dx[2]^4 -
6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2])*(sin(ξ)*dξ[1]^4 - 6*cos(ξ)*dξ[1]^2*dξ[2] - 3*sin(ξ)*dξ[2]^2 -
4*sin(ξ)*dξ[1]*dξ[3] + cos(ξ)*dξ[4]) + 56*(-(cos(ξ)*dξ[1]^3) - 3*sin(ξ)*dξ[1]*dξ[2] + cos(ξ)*dξ[3])*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] -
15*cos(x[2])*dx[2]*d2x[2]^2 - 10*cos(x[2])*dx[2]^2*d3x[2] - 10*sin(x[2])*d2x[2]*d3x[2] - 5*sin(x[2])*dx[2]*d4x[2] + cos(x[2])*d5x[2]) + 56*(-(cos(x[2])*dx[2]^3) -
3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2])*(cos(ξ)*dξ[1]^5 + 10*sin(ξ)*dξ[1]^3*dξ[2] - 15*cos(ξ)*dξ[1]*dξ[2]^2 - 10*cos(ξ)*dξ[1]^2*dξ[3] - 10*sin(ξ)*dξ[2]*dξ[3] -
5*sin(ξ)*dξ[1]*dξ[4] + cos(ξ)*dξ[5]) + 28*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2])*(-(sin(x[2])*dx[2]^6) + 15*cos(x[2])*dx[2]^4*d2x[2] + 45*sin(x[2])*dx[2]^2*d2x[2]^2 -
15*cos(x[2])*d2x[2]^3 + 20*sin(x[2])*dx[2]^3*d3x[2] - 60*cos(x[2])*dx[2]*d2x[2]*d3x[2] - 10*sin(x[2])*d3x[2]^2 - 15*cos(x[2])*dx[2]^2*d4x[2] - 15*sin(x[2])*d2x[2]*d4x[2] -
6*sin(x[2])*dx[2]*d5x[2] + cos(x[2])*d6x[2]) + 28*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(-(sin(ξ)*dξ[1]^6) + 15*cos(ξ)*dξ[1]^4*dξ[2] + 45*sin(ξ)*dξ[1]^2*dξ[2]^2 -
15*cos(ξ)*dξ[2]^3 + 20*sin(ξ)*dξ[1]^3*dξ[3] - 60*cos(ξ)*dξ[1]*dξ[2]*dξ[3] - 10*sin(ξ)*dξ[3]^2 - 15*cos(ξ)*dξ[1]^2*dξ[4] - 15*sin(ξ)*dξ[2]*dξ[4] - 6*sin(ξ)*dξ[1]*dξ[5] +
cos(ξ)*dξ[6]) + 8*cos(ξ)*dξ[1]*(-(cos(x[2])*dx[2]^7) - 21*sin(x[2])*dx[2]^5*d2x[2] + 105*cos(x[2])*dx[2]^3*d2x[2]^2 + 105*sin(x[2])*dx[2]*d2x[2]^3 +
35*cos(x[2])*dx[2]^4*d3x[2] + 210*sin(x[2])*dx[2]^2*d2x[2]*d3x[2] - 105*cos(x[2])*d2x[2]^2*d3x[2] - 70*cos(x[2])*dx[2]*d3x[2]^2 + 35*sin(x[2])*dx[2]^3*d4x[2] -
105*cos(x[2])*dx[2]*d2x[2]*d4x[2] - 35*sin(x[2])*d3x[2]*d4x[2] - 21*cos(x[2])*dx[2]^2*d5x[2] - 21*sin(x[2])*d2x[2]*d5x[2] - 7*sin(x[2])*dx[2]*d6x[2] + cos(x[2])*d7x[2]) +
8*cos(x[2])*dx[2]*(-(cos(ξ)*dξ[1]^7) - 21*sin(ξ)*dξ[1]^5*dξ[2] + 105*cos(ξ)*dξ[1]^3*dξ[2]^2 + 105*sin(ξ)*dξ[1]*dξ[2]^3 + 35*cos(ξ)*dξ[1]^4*dξ[3] +
210*sin(ξ)*dξ[1]^2*dξ[2]*dξ[3] - 105*cos(ξ)*dξ[2]^2*dξ[3] - 70*cos(ξ)*dξ[1]*dξ[3]^2 + 35*sin(ξ)*dξ[1]^3*dξ[4] - 105*cos(ξ)*dξ[1]*dξ[2]*dξ[4] - 35*sin(ξ)*dξ[3]*dξ[4] -
21*cos(ξ)*dξ[1]^2*dξ[5] - 21*sin(ξ)*dξ[2]*dξ[5] - 7*sin(ξ)*dξ[1]*dξ[6] + cos(ξ)*dξ[7]) + sin(ξ)*(sin(x[2])*dx[2]^8 - 28*cos(x[2])*dx[2]^6*d2x[2] -
210*sin(x[2])*dx[2]^4*d2x[2]^2 + 420*cos(x[2])*dx[2]^2*d2x[2]^3 + 105*sin(x[2])*d2x[2]^4 - 56*sin(x[2])*dx[2]^5*d3x[2] + 560*cos(x[2])*dx[2]^3*d2x[2]*d3x[2] +
840*sin(x[2])*dx[2]*d2x[2]^2*d3x[2] + 280*sin(x[2])*dx[2]^2*d3x[2]^2 - 280*cos(x[2])*d2x[2]*d3x[2]^2 + 70*cos(x[2])*dx[2]^4*d4x[2] + 420*sin(x[2])*dx[2]^2*d2x[2]*d4x[2] -
210*cos(x[2])*d2x[2]^2*d4x[2] - 280*cos(x[2])*dx[2]*d3x[2]*d4x[2] - 35*sin(x[2])*d4x[2]^2 + 56*sin(x[2])*dx[2]^3*d5x[2] - 168*cos(x[2])*dx[2]*d2x[2]*d5x[2] -
56*sin(x[2])*d3x[2]*d5x[2] - 28*cos(x[2])*dx[2]^2*d6x[2] - 28*sin(x[2])*d2x[2]*d6x[2] - 8*sin(x[2])*dx[2]*d7x[2] + cos(x[2])*d8x[2]) + sin(x[2])*(sin(ξ)*dξ[1]^8 -
28*cos(ξ)*dξ[1]^6*dξ[2] - 210*sin(ξ)*dξ[1]^4*dξ[2]^2 + 420*cos(ξ)*dξ[1]^2*dξ[2]^3 + 105*sin(ξ)*dξ[2]^4 - 56*sin(ξ)*dξ[1]^5*dξ[3] + 560*cos(ξ)*dξ[1]^3*dξ[2]*dξ[3] +
840*sin(ξ)*dξ[1]*dξ[2]^2*dξ[3] + 280*sin(ξ)*dξ[1]^2*dξ[3]^2 - 280*cos(ξ)*dξ[2]*dξ[3]^2 + 70*cos(ξ)*dξ[1]^4*dξ[4] + 420*sin(ξ)*dξ[1]^2*dξ[2]*dξ[4] -
210*cos(ξ)*dξ[2]^2*dξ[4] - 280*cos(ξ)*dξ[1]*dξ[3]*dξ[4] - 35*sin(ξ)*dξ[4]^2 + 56*sin(ξ)*dξ[1]^3*dξ[5] - 168*cos(ξ)*dξ[1]*dξ[2]*dξ[5] - 56*sin(ξ)*dξ[3]*dξ[5] -
28*cos(ξ)*dξ[1]^2*dξ[6] - 28*sin(ξ)*dξ[2]*dξ[6] - 8*sin(ξ)*dξ[1]*dξ[7] + cos(ξ)*dξ[8])) + sin(x[2])*sin(ξ)*(280*d2x[1]*d3x[1]^2*dρ[3] + 210*d2x[1]^2*dρ[3]*d4x[1] +
280*dx[1]*d3x[1]*dρ[3]*d4x[1] + 35*dρ[2]*d4x[1]^2 + 105*d2x[1]^4*dρ[4] + 840*dx[1]*d2x[1]^2*d3x[1]*dρ[4] + 280*dx[1]^2*d3x[1]^2*dρ[4] + 420*dx[1]^2*d2x[1]*d4x[1]*dρ[4] +
56*dρ[2]*d3x[1]*d5x[1] + 168*dx[1]*d2x[1]*dρ[3]*d5x[1] + 56*dx[1]^3*dρ[4]*d5x[1] + 420*dx[1]^2*d2x[1]^3*dρ[5] + 560*dx[1]^3*d2x[1]*d3x[1]*dρ[5] + 70*dx[1]^4*d4x[1]*dρ[5] +
28*d2x[1]*dρ[2]*d6x[1] + 28*dx[1]^2*dρ[3]*d6x[1] + 210*dx[1]^4*d2x[1]^2*dρ[6] + 56*dx[1]^5*d3x[1]*dρ[6] + 8*dx[1]*dρ[2]*d7x[1] + 28*dx[1]^6*d2x[1]*dρ[7] + dρ[1]*d8x[1] +
dx[1]^8*dρ[8])

zH(x::AbstractVector{Float64}, r::Float64) = (r - M) * cos(x[2])

dzH1(x::AbstractVector{Float64}, dx::AbstractVector{Float64}) = cos(x[2])*dx[1] - (-M + x[1])*sin(x[2])*dx[2]

dzH2(x::AbstractVector{Float64}, dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}) = -2*sin(x[2])*dx[1]*dx[2] + cos(x[2])*d2x[1] + (-M + x[1])*(-(cos(x[2])*dx[2]^2) - sin(x[2])*d2x[2])

dzH3(x::AbstractVector{Float64}, dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}) = -3*sin(x[2])*dx[2]*d2x[1] + 3*dx[1]*(-(cos(x[2])*dx[2]^2) -
sin(x[2])*d2x[2]) + cos(x[2])*d3x[1] + (-M + x[1])*(sin(x[2])*dx[2]^3 - 3*cos(x[2])*dx[2]*d2x[2] - sin(x[2])*d3x[2])

dzH4(x::AbstractVector{Float64}, dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, d4x::AbstractVector{Float64}) = 6*d2x[1]*(-(cos(x[2])*dx[2]^2) - sin(x[2])*d2x[2]) -
4*sin(x[2])*dx[2]*d3x[1] + 4*dx[1]*(sin(x[2])*dx[2]^3 - 3*cos(x[2])*dx[2]*d2x[2] - sin(x[2])*d3x[2]) + cos(x[2])*d4x[1] + (-M + x[1])*(cos(x[2])*dx[2]^4 +
6*sin(x[2])*dx[2]^2*d2x[2] - 3*cos(x[2])*d2x[2]^2 - 4*cos(x[2])*dx[2]*d3x[2] - sin(x[2])*d4x[2])

dzH5(x::AbstractVector{Float64}, dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, d4x::AbstractVector{Float64}, d5x::AbstractVector{Float64}) = 10*(-(cos(x[2])*dx[2]^2) -
sin(x[2])*d2x[2])*d3x[1] + 10*d2x[1]*(sin(x[2])*dx[2]^3 - 3*cos(x[2])*dx[2]*d2x[2] - sin(x[2])*d3x[2]) - 5*sin(x[2])*dx[2]*d4x[1] + 5*dx[1]*(cos(x[2])*dx[2]^4 +
6*sin(x[2])*dx[2]^2*d2x[2] - 3*cos(x[2])*d2x[2]^2 - 4*cos(x[2])*dx[2]*d3x[2] - sin(x[2])*d4x[2]) + cos(x[2])*d5x[1] + (-M + x[1])*(-(sin(x[2])*dx[2]^5) +
10*cos(x[2])*dx[2]^3*d2x[2] + 15*sin(x[2])*dx[2]*d2x[2]^2 + 10*sin(x[2])*dx[2]^2*d3x[2] - 10*cos(x[2])*d2x[2]*d3x[2] - 5*cos(x[2])*dx[2]*d4x[2] - sin(x[2])*d5x[2])

dzH6(x::AbstractVector{Float64}, dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, d4x::AbstractVector{Float64}, d5x::AbstractVector{Float64}, d6x::AbstractVector{Float64}) = 20*d3x[1]*
(sin(x[2])*dx[2]^3 - 3*cos(x[2])*dx[2]*d2x[2] - sin(x[2])*d3x[2]) + 15*(-(cos(x[2])*dx[2]^2) - sin(x[2])*d2x[2])*d4x[1] + 15*d2x[1]*(cos(x[2])*dx[2]^4 +
6*sin(x[2])*dx[2]^2*d2x[2] - 3*cos(x[2])*d2x[2]^2 - 4*cos(x[2])*dx[2]*d3x[2] - sin(x[2])*d4x[2]) - 6*sin(x[2])*dx[2]*d5x[1] + 6*dx[1]*(-(sin(x[2])*dx[2]^5) +
10*cos(x[2])*dx[2]^3*d2x[2] + 15*sin(x[2])*dx[2]*d2x[2]^2 + 10*sin(x[2])*dx[2]^2*d3x[2] - 10*cos(x[2])*d2x[2]*d3x[2] - 5*cos(x[2])*dx[2]*d4x[2] - sin(x[2])*d5x[2]) +
cos(x[2])*d6x[1] + (-M + x[1])*(-(cos(x[2])*dx[2]^6) - 15*sin(x[2])*dx[2]^4*d2x[2] + 45*cos(x[2])*dx[2]^2*d2x[2]^2 + 15*sin(x[2])*d2x[2]^3 + 20*cos(x[2])*dx[2]^3*d3x[2] +
60*sin(x[2])*dx[2]*d2x[2]*d3x[2] - 10*cos(x[2])*d3x[2]^2 + 15*sin(x[2])*dx[2]^2*d4x[2] - 15*cos(x[2])*d2x[2]*d4x[2] - 6*cos(x[2])*dx[2]*d5x[2] - sin(x[2])*d6x[2])

dzH7(x::AbstractVector{Float64}, dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, d4x::AbstractVector{Float64}, d5x::AbstractVector{Float64}, d6x::AbstractVector{Float64}, d7x::AbstractVector{Float64}) = 35*(sin(x[2])*dx[2]^3 - 3*cos(x[2])*dx[2]*d2x[2] - sin(x[2])*d3x[2])*d4x[1] + 35*d3x[1]*(cos(x[2])*dx[2]^4 + 6*sin(x[2])*dx[2]^2*d2x[2] - 3*cos(x[2])*d2x[2]^2 -
4*cos(x[2])*dx[2]*d3x[2] - sin(x[2])*d4x[2]) + 21*(-(cos(x[2])*dx[2]^2) - sin(x[2])*d2x[2])*d5x[1] + 21*d2x[1]*(-(sin(x[2])*dx[2]^5) + 10*cos(x[2])*dx[2]^3*d2x[2] +
15*sin(x[2])*dx[2]*d2x[2]^2 + 10*sin(x[2])*dx[2]^2*d3x[2] - 10*cos(x[2])*d2x[2]*d3x[2] - 5*cos(x[2])*dx[2]*d4x[2] - sin(x[2])*d5x[2]) - 7*sin(x[2])*dx[2]*d6x[1] +
7*dx[1]*(-(cos(x[2])*dx[2]^6) - 15*sin(x[2])*dx[2]^4*d2x[2] + 45*cos(x[2])*dx[2]^2*d2x[2]^2 + 15*sin(x[2])*d2x[2]^3 + 20*cos(x[2])*dx[2]^3*d3x[2] +
60*sin(x[2])*dx[2]*d2x[2]*d3x[2] - 10*cos(x[2])*d3x[2]^2 + 15*sin(x[2])*dx[2]^2*d4x[2] - 15*cos(x[2])*d2x[2]*d4x[2] - 6*cos(x[2])*dx[2]*d5x[2] - sin(x[2])*d6x[2]) +
cos(x[2])*d7x[1] + (-M + x[1])*(sin(x[2])*dx[2]^7 - 21*cos(x[2])*dx[2]^5*d2x[2] - 105*sin(x[2])*dx[2]^3*d2x[2]^2 + 105*cos(x[2])*dx[2]*d2x[2]^3 -
35*sin(x[2])*dx[2]^4*d3x[2] + 210*cos(x[2])*dx[2]^2*d2x[2]*d3x[2] + 105*sin(x[2])*d2x[2]^2*d3x[2] + 70*sin(x[2])*dx[2]*d3x[2]^2 + 35*cos(x[2])*dx[2]^3*d4x[2] +
105*sin(x[2])*dx[2]*d2x[2]*d4x[2] - 35*cos(x[2])*d3x[2]*d4x[2] + 21*sin(x[2])*dx[2]^2*d5x[2] - 21*cos(x[2])*d2x[2]*d5x[2] - 7*cos(x[2])*dx[2]*d6x[2] - sin(x[2])*d7x[2])

dzH8(x::AbstractVector{Float64}, dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, d4x::AbstractVector{Float64}, d5x::AbstractVector{Float64}, d6x::AbstractVector{Float64}, d7x::AbstractVector{Float64}, d8x::AbstractVector{Float64}) = 70*d4x[1]*(cos(x[2])*dx[2]^4 + 6*sin(x[2])*dx[2]^2*d2x[2] - 3*cos(x[2])*d2x[2]^2 - 4*cos(x[2])*dx[2]*d3x[2] - sin(x[2])*d4x[2]) +
56*(sin(x[2])*dx[2]^3 - 3*cos(x[2])*dx[2]*d2x[2] - sin(x[2])*d3x[2])*d5x[1] + 56*d3x[1]*(-(sin(x[2])*dx[2]^5) + 10*cos(x[2])*dx[2]^3*d2x[2] + 15*sin(x[2])*dx[2]*d2x[2]^2 +
10*sin(x[2])*dx[2]^2*d3x[2] - 10*cos(x[2])*d2x[2]*d3x[2] - 5*cos(x[2])*dx[2]*d4x[2] - sin(x[2])*d5x[2]) + 28*(-(cos(x[2])*dx[2]^2) - sin(x[2])*d2x[2])*d6x[1] +
28*d2x[1]*(-(cos(x[2])*dx[2]^6) - 15*sin(x[2])*dx[2]^4*d2x[2] + 45*cos(x[2])*dx[2]^2*d2x[2]^2 + 15*sin(x[2])*d2x[2]^3 + 20*cos(x[2])*dx[2]^3*d3x[2] +
60*sin(x[2])*dx[2]*d2x[2]*d3x[2] - 10*cos(x[2])*d3x[2]^2 + 15*sin(x[2])*dx[2]^2*d4x[2] - 15*cos(x[2])*d2x[2]*d4x[2] - 6*cos(x[2])*dx[2]*d5x[2] - sin(x[2])*d6x[2]) -
8*sin(x[2])*dx[2]*d7x[1] + 8*dx[1]*(sin(x[2])*dx[2]^7 - 21*cos(x[2])*dx[2]^5*d2x[2] - 105*sin(x[2])*dx[2]^3*d2x[2]^2 + 105*cos(x[2])*dx[2]*d2x[2]^3 -
35*sin(x[2])*dx[2]^4*d3x[2] + 210*cos(x[2])*dx[2]^2*d2x[2]*d3x[2] + 105*sin(x[2])*d2x[2]^2*d3x[2] + 70*sin(x[2])*dx[2]*d3x[2]^2 + 35*cos(x[2])*dx[2]^3*d4x[2] +
105*sin(x[2])*dx[2]*d2x[2]*d4x[2] - 35*cos(x[2])*d3x[2]*d4x[2] + 21*sin(x[2])*dx[2]^2*d5x[2] - 21*cos(x[2])*d2x[2]*d5x[2] - 7*cos(x[2])*dx[2]*d6x[2] - sin(x[2])*d7x[2]) +
cos(x[2])*d8x[1] + (-M + x[1])*(cos(x[2])*dx[2]^8 + 28*sin(x[2])*dx[2]^6*d2x[2] - 210*cos(x[2])*dx[2]^4*d2x[2]^2 - 420*sin(x[2])*dx[2]^2*d2x[2]^3 + 105*cos(x[2])*d2x[2]^4 -
56*cos(x[2])*dx[2]^5*d3x[2] - 560*sin(x[2])*dx[2]^3*d2x[2]*d3x[2] + 840*cos(x[2])*dx[2]*d2x[2]^2*d3x[2] + 280*cos(x[2])*dx[2]^2*d3x[2]^2 + 280*sin(x[2])*d2x[2]*d3x[2]^2 -
70*sin(x[2])*dx[2]^4*d4x[2] + 420*cos(x[2])*dx[2]^2*d2x[2]*d4x[2] + 210*sin(x[2])*d2x[2]^2*d4x[2] + 280*sin(x[2])*dx[2]*d3x[2]*d4x[2] - 35*cos(x[2])*d4x[2]^2 +
56*cos(x[2])*dx[2]^3*d5x[2] + 168*sin(x[2])*dx[2]*d2x[2]*d5x[2] - 56*cos(x[2])*d3x[2]*d5x[2] + 28*sin(x[2])*dx[2]^2*d6x[2] - 28*cos(x[2])*d2x[2]*d6x[2] -
8*cos(x[2])*dx[2]*d7x[2] - sin(x[2])*d8x[2])



dxH9(x::AbstractVector{Float64}, dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, d4x::AbstractVector{Float64}, d5x::AbstractVector{Float64}, d6x::AbstractVector{Float64}, d7x::AbstractVector{Float64}, d8x::AbstractVector{Float64}, d9x::AbstractVector{Float64}, ξ::Float64, dξ::AbstractVector{Float64}, ρ::Float64, dρ::AbstractVector{Float64}) = 126*(3*d2x[1]^2*dρ[2] + 4*dx[1]*dρ[2]*d3x[1] + 6*dx[1]^2*d2x[1]*dρ[3] + dρ[1]*d4x[1] + dx[1]^4*dρ[4])*(10*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2])*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) + 10*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(sin(ξ)*dξ[1]^3 - 3*cos(ξ)*dξ[1]*dξ[2] - sin(ξ)*dξ[3]) - 5*sin(ξ)*dξ[1]*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + 5*cos(x[2])*dx[2]*(cos(ξ)*dξ[1]^4 + 6*sin(ξ)*dξ[1]^2*dξ[2] - 3*cos(ξ)*dξ[2]^2 - 4*cos(ξ)*dξ[1]*dξ[3] - sin(ξ)*dξ[4]) + cos(ξ)*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] - 15*cos(x[2])*dx[2]*d2x[2]^2 - 10*cos(x[2])*dx[2]^2*d3x[2] - 10*sin(x[2])*d2x[2]*d3x[2] - 5*sin(x[2])*dx[2]*d4x[2] + cos(x[2])*d5x[2]) + sin(x[2])*(-(sin(ξ)*dξ[1]^5) + 10*cos(ξ)*dξ[1]^3*dξ[2] + 15*sin(ξ)*dξ[1]*dξ[2]^2 + 10*sin(ξ)*dξ[1]^2*dξ[3] - 10*cos(ξ)*dξ[2]*dξ[3] - 5*cos(ξ)*dξ[1]*dξ[4] - sin(ξ)*dξ[5])) + 126*(6*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2]) - 4*sin(ξ)*dξ[1]*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) + 4*cos(x[2])*dx[2]*(sin(ξ)*dξ[1]^3 - 3*cos(ξ)*dξ[1]*dξ[2] - sin(ξ)*dξ[3]) + cos(ξ)*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + sin(x[2])*(cos(ξ)*dξ[1]^4 + 6*sin(ξ)*dξ[1]^2*dξ[2] - 3*cos(ξ)*dξ[2]^2 - 4*cos(ξ)*dξ[1]*dξ[3] - sin(ξ)*dξ[4]))*(10*d2x[1]*dρ[2]*d3x[1] + 15*dx[1]*d2x[1]^2*dρ[3] + 10*dx[1]^2*d3x[1]*dρ[3] + 5*dx[1]*dρ[2]*d4x[1] + 10*dx[1]^3*d2x[1]*dρ[4] + dρ[1]*d5x[1] + dx[1]^5*dρ[5]) + 84*(3*dx[1]*d2x[1]*dρ[2] + dρ[1]*d3x[1] + dx[1]^3*dρ[3])*(20*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2])*(sin(ξ)*dξ[1]^3 - 3*cos(ξ)*dξ[1]*dξ[2] - sin(ξ)*dξ[3]) + 15*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2])*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + 15*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(cos(ξ)*dξ[1]^4 + 6*sin(ξ)*dξ[1]^2*dξ[2] - 3*cos(ξ)*dξ[2]^2 - 4*cos(ξ)*dξ[1]*dξ[3] - sin(ξ)*dξ[4]) - 6*sin(ξ)*dξ[1]*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] - 15*cos(x[2])*dx[2]*d2x[2]^2 - 10*cos(x[2])*dx[2]^2*d3x[2] - 10*sin(x[2])*d2x[2]*d3x[2] - 5*sin(x[2])*dx[2]*d4x[2] + cos(x[2])*d5x[2]) + 6*cos(x[2])*dx[2]*(-(sin(ξ)*dξ[1]^5) + 10*cos(ξ)*dξ[1]^3*dξ[2] + 15*sin(ξ)*dξ[1]*dξ[2]^2 + 10*sin(ξ)*dξ[1]^2*dξ[3] - 10*cos(ξ)*dξ[2]*dξ[3] - 5*cos(ξ)*dξ[1]*dξ[4] - sin(ξ)*dξ[5]) + cos(ξ)*(-(sin(x[2])*dx[2]^6) + 15*cos(x[2])*dx[2]^4*d2x[2] + 45*sin(x[2])*dx[2]^2*d2x[2]^2 - 15*cos(x[2])*d2x[2]^3 + 20*sin(x[2])*dx[2]^3*d3x[2] - 60*cos(x[2])*dx[2]*d2x[2]*d3x[2] - 10*sin(x[2])*d3x[2]^2 - 15*cos(x[2])*dx[2]^2*d4x[2] - 15*sin(x[2])*d2x[2]*d4x[2] - 6*sin(x[2])*dx[2]*d5x[2] + cos(x[2])*d6x[2]) + sin(x[2])*(-(cos(ξ)*dξ[1]^6) - 15*sin(ξ)*dξ[1]^4*dξ[2] + 45*cos(ξ)*dξ[1]^2*dξ[2]^2 + 15*sin(ξ)*dξ[2]^3 + 20*cos(ξ)*dξ[1]^3*dξ[3] + 60*sin(ξ)*dξ[1]*dξ[2]*dξ[3] - 10*cos(ξ)*dξ[3]^2 + 15*sin(ξ)*dξ[1]^2*dξ[4] - 15*cos(ξ)*dξ[2]*dξ[4] - 6*cos(ξ)*dξ[1]*dξ[5] - sin(ξ)*dξ[6])) + 84*(-3*sin(ξ)*dξ[1]*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) + 3*cos(x[2])*dx[2]*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2]) + cos(ξ)*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) + sin(x[2])*(sin(ξ)*dξ[1]^3 - 3*cos(ξ)*dξ[1]*dξ[2] - sin(ξ)*dξ[3]))*(10*dρ[2]*d3x[1]^2 + 15*d2x[1]^3*dρ[3] + 60*dx[1]*d2x[1]*d3x[1]*dρ[3] + 15*d2x[1]*dρ[2]*d4x[1] + 15*dx[1]^2*dρ[3]*d4x[1] + 45*dx[1]^2*d2x[1]^2*dρ[4] + 20*dx[1]^3*d3x[1]*dρ[4] + 6*dx[1]*dρ[2]*d5x[1] + 15*dx[1]^4*d2x[1]*dρ[5] + dρ[1]*d6x[1] + dx[1]^6*dρ[6]) + 36*(dρ[1]*d2x[1] + dx[1]^2*dρ[2])*(35*(sin(ξ)*dξ[1]^3 - 3*cos(ξ)*dξ[1]*dξ[2] - sin(ξ)*dξ[3])*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + 35*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2])*(cos(ξ)*dξ[1]^4 + 6*sin(ξ)*dξ[1]^2*dξ[2] - 3*cos(ξ)*dξ[2]^2 - 4*cos(ξ)*dξ[1]*dξ[3] - sin(ξ)*dξ[4]) + 21*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2])*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] - 15*cos(x[2])*dx[2]*d2x[2]^2 - 10*cos(x[2])*dx[2]^2*d3x[2] - 10*sin(x[2])*d2x[2]*d3x[2] - 5*sin(x[2])*dx[2]*d4x[2] + cos(x[2])*d5x[2]) + 21*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(-(sin(ξ)*dξ[1]^5) + 10*cos(ξ)*dξ[1]^3*dξ[2] + 15*sin(ξ)*dξ[1]*dξ[2]^2 + 10*sin(ξ)*dξ[1]^2*dξ[3] - 10*cos(ξ)*dξ[2]*dξ[3] - 5*cos(ξ)*dξ[1]*dξ[4] - sin(ξ)*dξ[5]) - 7*sin(ξ)*dξ[1]*(-(sin(x[2])*dx[2]^6) + 15*cos(x[2])*dx[2]^4*d2x[2] + 45*sin(x[2])*dx[2]^2*d2x[2]^2 - 15*cos(x[2])*d2x[2]^3 + 20*sin(x[2])*dx[2]^3*d3x[2] - 60*cos(x[2])*dx[2]*d2x[2]*d3x[2] - 10*sin(x[2])*d3x[2]^2 - 15*cos(x[2])*dx[2]^2*d4x[2] - 15*sin(x[2])*d2x[2]*d4x[2] - 6*sin(x[2])*dx[2]*d5x[2] + cos(x[2])*d6x[2]) + 7*cos(x[2])*dx[2]*(-(cos(ξ)*dξ[1]^6) - 15*sin(ξ)*dξ[1]^4*dξ[2] + 45*cos(ξ)*dξ[1]^2*dξ[2]^2 + 15*sin(ξ)*dξ[2]^3 + 20*cos(ξ)*dξ[1]^3*dξ[3] + 60*sin(ξ)*dξ[1]*dξ[2]*dξ[3] - 10*cos(ξ)*dξ[3]^2 + 15*sin(ξ)*dξ[1]^2*dξ[4] - 15*cos(ξ)*dξ[2]*dξ[4] - 6*cos(ξ)*dξ[1]*dξ[5] - sin(ξ)*dξ[6]) + cos(ξ)*(-(cos(x[2])*dx[2]^7) - 21*sin(x[2])*dx[2]^5*d2x[2] + 105*cos(x[2])*dx[2]^3*d2x[2]^2 + 105*sin(x[2])*dx[2]*d2x[2]^3 + 35*cos(x[2])*dx[2]^4*d3x[2] + 210*sin(x[2])*dx[2]^2*d2x[2]*d3x[2] - 105*cos(x[2])*d2x[2]^2*d3x[2] - 70*cos(x[2])*dx[2]*d3x[2]^2 + 35*sin(x[2])*dx[2]^3*d4x[2] - 105*cos(x[2])*dx[2]*d2x[2]*d4x[2] - 35*sin(x[2])*d3x[2]*d4x[2] - 21*cos(x[2])*dx[2]^2*d5x[2] - 21*sin(x[2])*d2x[2]*d5x[2] - 7*sin(x[2])*dx[2]*d6x[2] + cos(x[2])*d7x[2]) + sin(x[2])*(sin(ξ)*dξ[1]^7 - 21*cos(ξ)*dξ[1]^5*dξ[2] - 105*sin(ξ)*dξ[1]^3*dξ[2]^2 + 105*cos(ξ)*dξ[1]*dξ[2]^3 - 35*sin(ξ)*dξ[1]^4*dξ[3] + 210*cos(ξ)*dξ[1]^2*dξ[2]*dξ[3] + 105*sin(ξ)*dξ[2]^2*dξ[3] + 70*sin(ξ)*dξ[1]*dξ[3]^2 + 35*cos(ξ)*dξ[1]^3*dξ[4] + 105*sin(ξ)*dξ[1]*dξ[2]*dξ[4] - 35*cos(ξ)*dξ[3]*dξ[4] + 21*sin(ξ)*dξ[1]^2*dξ[5] - 21*cos(ξ)*dξ[2]*dξ[5] - 7*cos(ξ)*dξ[1]*dξ[6] - sin(ξ)*dξ[7])) + 36*(-2*cos(x[2])*sin(ξ)*dx[2]*dξ[1] + cos(ξ)*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) + sin(x[2])*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2]))*(105*d2x[1]^2*d3x[1]*dρ[3] + 70*dx[1]*d3x[1]^2*dρ[3] + 35*dρ[2]*d3x[1]*d4x[1] + 105*dx[1]*d2x[1]*dρ[3]*d4x[1] + 105*dx[1]*d2x[1]^3*dρ[4] + 210*dx[1]^2*d2x[1]*d3x[1]*dρ[4] + 35*dx[1]^3*d4x[1]*dρ[4] + 21*d2x[1]*dρ[2]*d5x[1] + 21*dx[1]^2*dρ[3]*d5x[1] + 105*dx[1]^3*d2x[1]^2*dρ[5] + 35*dx[1]^4*d3x[1]*dρ[5] + 7*dx[1]*dρ[2]*d6x[1] + 21*dx[1]^5*d2x[1]*dρ[6] + dρ[1]*d7x[1] + dx[1]^7*dρ[7]) + 9*dx[1]*dρ[1]*(70*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2])*(cos(ξ)*dξ[1]^4 + 6*sin(ξ)*dξ[1]^2*dξ[2] - 3*cos(ξ)*dξ[2]^2 - 4*cos(ξ)*dξ[1]*dξ[3] - sin(ξ)*dξ[4]) + 56*(sin(ξ)*dξ[1]^3 - 3*cos(ξ)*dξ[1]*dξ[2] - sin(ξ)*dξ[3])*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] - 15*cos(x[2])*dx[2]*d2x[2]^2 - 10*cos(x[2])*dx[2]^2*d3x[2] - 10*sin(x[2])*d2x[2]*d3x[2] - 5*sin(x[2])*dx[2]*d4x[2] + cos(x[2])*d5x[2]) + 56*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2])*(-(sin(ξ)*dξ[1]^5) + 10*cos(ξ)*dξ[1]^3*dξ[2] + 15*sin(ξ)*dξ[1]*dξ[2]^2 + 10*sin(ξ)*dξ[1]^2*dξ[3] - 10*cos(ξ)*dξ[2]*dξ[3] - 5*cos(ξ)*dξ[1]*dξ[4] - sin(ξ)*dξ[5]) + 28*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2])*(-(sin(x[2])*dx[2]^6) + 15*cos(x[2])*dx[2]^4*d2x[2] + 45*sin(x[2])*dx[2]^2*d2x[2]^2 - 15*cos(x[2])*d2x[2]^3 + 20*sin(x[2])*dx[2]^3*d3x[2] - 60*cos(x[2])*dx[2]*d2x[2]*d3x[2] - 10*sin(x[2])*d3x[2]^2 - 15*cos(x[2])*dx[2]^2*d4x[2] - 15*sin(x[2])*d2x[2]*d4x[2] - 6*sin(x[2])*dx[2]*d5x[2] + cos(x[2])*d6x[2]) + 28*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(-(cos(ξ)*dξ[1]^6) - 15*sin(ξ)*dξ[1]^4*dξ[2] + 45*cos(ξ)*dξ[1]^2*dξ[2]^2 + 15*sin(ξ)*dξ[2]^3 + 20*cos(ξ)*dξ[1]^3*dξ[3] + 60*sin(ξ)*dξ[1]*dξ[2]*dξ[3] - 10*cos(ξ)*dξ[3]^2 + 15*sin(ξ)*dξ[1]^2*dξ[4] - 15*cos(ξ)*dξ[2]*dξ[4] - 6*cos(ξ)*dξ[1]*dξ[5] - sin(ξ)*dξ[6]) - 8*sin(ξ)*dξ[1]*(-(cos(x[2])*dx[2]^7) - 21*sin(x[2])*dx[2]^5*d2x[2] + 105*cos(x[2])*dx[2]^3*d2x[2]^2 + 105*sin(x[2])*dx[2]*d2x[2]^3 + 35*cos(x[2])*dx[2]^4*d3x[2] + 210*sin(x[2])*dx[2]^2*d2x[2]*d3x[2] - 105*cos(x[2])*d2x[2]^2*d3x[2] - 70*cos(x[2])*dx[2]*d3x[2]^2 + 35*sin(x[2])*dx[2]^3*d4x[2] - 105*cos(x[2])*dx[2]*d2x[2]*d4x[2] - 35*sin(x[2])*d3x[2]*d4x[2] - 21*cos(x[2])*dx[2]^2*d5x[2] - 21*sin(x[2])*d2x[2]*d5x[2] - 7*sin(x[2])*dx[2]*d6x[2] + cos(x[2])*d7x[2]) + 8*cos(x[2])*dx[2]*(sin(ξ)*dξ[1]^7 - 21*cos(ξ)*dξ[1]^5*dξ[2] - 105*sin(ξ)*dξ[1]^3*dξ[2]^2 + 105*cos(ξ)*dξ[1]*dξ[2]^3 - 35*sin(ξ)*dξ[1]^4*dξ[3] + 210*cos(ξ)*dξ[1]^2*dξ[2]*dξ[3] + 105*sin(ξ)*dξ[2]^2*dξ[3] + 70*sin(ξ)*dξ[1]*dξ[3]^2 + 35*cos(ξ)*dξ[1]^3*dξ[4] + 105*sin(ξ)*dξ[1]*dξ[2]*dξ[4] - 35*cos(ξ)*dξ[3]*dξ[4] + 21*sin(ξ)*dξ[1]^2*dξ[5] - 21*cos(ξ)*dξ[2]*dξ[5] - 7*cos(ξ)*dξ[1]*dξ[6] - sin(ξ)*dξ[7]) + cos(ξ)*(sin(x[2])*dx[2]^8 - 28*cos(x[2])*dx[2]^6*d2x[2] - 210*sin(x[2])*dx[2]^4*d2x[2]^2 + 420*cos(x[2])*dx[2]^2*d2x[2]^3 + 105*sin(x[2])*d2x[2]^4 - 56*sin(x[2])*dx[2]^5*d3x[2] + 560*cos(x[2])*dx[2]^3*d2x[2]*d3x[2] + 840*sin(x[2])*dx[2]*d2x[2]^2*d3x[2] + 280*sin(x[2])*dx[2]^2*d3x[2]^2 - 280*cos(x[2])*d2x[2]*d3x[2]^2 + 70*cos(x[2])*dx[2]^4*d4x[2] + 420*sin(x[2])*dx[2]^2*d2x[2]*d4x[2] - 210*cos(x[2])*d2x[2]^2*d4x[2] - 280*cos(x[2])*dx[2]*d3x[2]*d4x[2] - 35*sin(x[2])*d4x[2]^2 + 56*sin(x[2])*dx[2]^3*d5x[2] - 168*cos(x[2])*dx[2]*d2x[2]*d5x[2] - 56*sin(x[2])*d3x[2]*d5x[2] - 28*cos(x[2])*dx[2]^2*d6x[2] - 28*sin(x[2])*d2x[2]*d6x[2] - 8*sin(x[2])*dx[2]*d7x[2] + cos(x[2])*d8x[2]) + sin(x[2])*(cos(ξ)*dξ[1]^8 + 28*sin(ξ)*dξ[1]^6*dξ[2] - 210*cos(ξ)*dξ[1]^4*dξ[2]^2 - 420*sin(ξ)*dξ[1]^2*dξ[2]^3 + 105*cos(ξ)*dξ[2]^4 - 56*cos(ξ)*dξ[1]^5*dξ[3] - 560*sin(ξ)*dξ[1]^3*dξ[2]*dξ[3] + 840*cos(ξ)*dξ[1]*dξ[2]^2*dξ[3] + 280*cos(ξ)*dξ[1]^2*dξ[3]^2 + 280*sin(ξ)*dξ[2]*dξ[3]^2 - 70*sin(ξ)*dξ[1]^4*dξ[4] + 420*cos(ξ)*dξ[1]^2*dξ[2]*dξ[4] + 210*sin(ξ)*dξ[2]^2*dξ[4] + 280*sin(ξ)*dξ[1]*dξ[3]*dξ[4] - 35*cos(ξ)*dξ[4]^2 + 56*cos(ξ)*dξ[1]^3*dξ[5] + 168*sin(ξ)*dξ[1]*dξ[2]*dξ[5] - 56*cos(ξ)*dξ[3]*dξ[5] + 28*sin(ξ)*dξ[1]^2*dξ[6] - 28*cos(ξ)*dξ[2]*dξ[6] - 8*cos(ξ)*dξ[1]*dξ[7] - sin(ξ)*dξ[8])) + 9*(cos(x[2])*cos(ξ)*dx[2] - sin(x[2])*sin(ξ)*dξ[1])*(280*d2x[1]*d3x[1]^2*dρ[3] + 210*d2x[1]^2*dρ[3]*d4x[1] + 280*dx[1]*d3x[1]*dρ[3]*d4x[1] + 35*dρ[2]*d4x[1]^2 + 105*d2x[1]^4*dρ[4] + 840*dx[1]*d2x[1]^2*d3x[1]*dρ[4] + 280*dx[1]^2*d3x[1]^2*dρ[4] + 420*dx[1]^2*d2x[1]*d4x[1]*dρ[4] + 56*dρ[2]*d3x[1]*d5x[1] + 168*dx[1]*d2x[1]*dρ[3]*d5x[1] + 56*dx[1]^3*dρ[4]*d5x[1] + 420*dx[1]^2*d2x[1]^3*dρ[5] + 560*dx[1]^3*d2x[1]*d3x[1]*dρ[5] + 70*dx[1]^4*d4x[1]*dρ[5] + 28*d2x[1]*dρ[2]*d6x[1] + 28*dx[1]^2*dρ[3]*d6x[1] + 210*dx[1]^4*d2x[1]^2*dρ[6] + 56*dx[1]^5*d3x[1]*dρ[6] + 8*dx[1]*dρ[2]*d7x[1] + 28*dx[1]^6*d2x[1]*dρ[7] + dρ[1]*d8x[1] + dx[1]^8*dρ[8]) + ρ*(126*(cos(ξ)*dξ[1]^4 + 6*sin(ξ)*dξ[1]^2*dξ[2] - 3*cos(ξ)*dξ[2]^2 - 4*cos(ξ)*dξ[1]*dξ[3] - sin(ξ)*dξ[4])*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] - 15*cos(x[2])*dx[2]*d2x[2]^2 - 10*cos(x[2])*dx[2]^2*d3x[2] - 10*sin(x[2])*d2x[2]*d3x[2] - 5*sin(x[2])*dx[2]*d4x[2] + cos(x[2])*d5x[2]) + 126*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2])*(-(sin(ξ)*dξ[1]^5) + 10*cos(ξ)*dξ[1]^3*dξ[2] + 15*sin(ξ)*dξ[1]*dξ[2]^2 + 10*sin(ξ)*dξ[1]^2*dξ[3] - 10*cos(ξ)*dξ[2]*dξ[3] - 5*cos(ξ)*dξ[1]*dξ[4] - sin(ξ)*dξ[5]) + 84*(sin(ξ)*dξ[1]^3 - 3*cos(ξ)*dξ[1]*dξ[2] - sin(ξ)*dξ[3])*(-(sin(x[2])*dx[2]^6) + 15*cos(x[2])*dx[2]^4*d2x[2] + 45*sin(x[2])*dx[2]^2*d2x[2]^2 - 15*cos(x[2])*d2x[2]^3 + 20*sin(x[2])*dx[2]^3*d3x[2] - 60*cos(x[2])*dx[2]*d2x[2]*d3x[2] - 10*sin(x[2])*d3x[2]^2 - 15*cos(x[2])*dx[2]^2*d4x[2] - 15*sin(x[2])*d2x[2]*d4x[2] - 6*sin(x[2])*dx[2]*d5x[2] + cos(x[2])*d6x[2]) + 84*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2])*(-(cos(ξ)*dξ[1]^6) - 15*sin(ξ)*dξ[1]^4*dξ[2] + 45*cos(ξ)*dξ[1]^2*dξ[2]^2 + 15*sin(ξ)*dξ[2]^3 + 20*cos(ξ)*dξ[1]^3*dξ[3] + 60*sin(ξ)*dξ[1]*dξ[2]*dξ[3] - 10*cos(ξ)*dξ[3]^2 + 15*sin(ξ)*dξ[1]^2*dξ[4] - 15*cos(ξ)*dξ[2]*dξ[4] - 6*cos(ξ)*dξ[1]*dξ[5] - sin(ξ)*dξ[6]) + 36*(-(cos(ξ)*dξ[1]^2) - sin(ξ)*dξ[2])*(-(cos(x[2])*dx[2]^7) - 21*sin(x[2])*dx[2]^5*d2x[2] + 105*cos(x[2])*dx[2]^3*d2x[2]^2 + 105*sin(x[2])*dx[2]*d2x[2]^3 + 35*cos(x[2])*dx[2]^4*d3x[2] + 210*sin(x[2])*dx[2]^2*d2x[2]*d3x[2] - 105*cos(x[2])*d2x[2]^2*d3x[2] - 70*cos(x[2])*dx[2]*d3x[2]^2 + 35*sin(x[2])*dx[2]^3*d4x[2] - 105*cos(x[2])*dx[2]*d2x[2]*d4x[2] - 35*sin(x[2])*d3x[2]*d4x[2] - 21*cos(x[2])*dx[2]^2*d5x[2] - 21*sin(x[2])*d2x[2]*d5x[2] - 7*sin(x[2])*dx[2]*d6x[2] + cos(x[2])*d7x[2]) + 36*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(sin(ξ)*dξ[1]^7 - 21*cos(ξ)*dξ[1]^5*dξ[2] - 105*sin(ξ)*dξ[1]^3*dξ[2]^2 + 105*cos(ξ)*dξ[1]*dξ[2]^3 - 35*sin(ξ)*dξ[1]^4*dξ[3] + 210*cos(ξ)*dξ[1]^2*dξ[2]*dξ[3] + 105*sin(ξ)*dξ[2]^2*dξ[3] + 70*sin(ξ)*dξ[1]*dξ[3]^2 + 35*cos(ξ)*dξ[1]^3*dξ[4] + 105*sin(ξ)*dξ[1]*dξ[2]*dξ[4] - 35*cos(ξ)*dξ[3]*dξ[4] + 21*sin(ξ)*dξ[1]^2*dξ[5] - 21*cos(ξ)*dξ[2]*dξ[5] - 7*cos(ξ)*dξ[1]*dξ[6] - sin(ξ)*dξ[7]) - 9*sin(ξ)*dξ[1]*(sin(x[2])*dx[2]^8 - 28*cos(x[2])*dx[2]^6*d2x[2] - 210*sin(x[2])*dx[2]^4*d2x[2]^2 + 420*cos(x[2])*dx[2]^2*d2x[2]^3 + 105*sin(x[2])*d2x[2]^4 - 56*sin(x[2])*dx[2]^5*d3x[2] + 560*cos(x[2])*dx[2]^3*d2x[2]*d3x[2] + 840*sin(x[2])*dx[2]*d2x[2]^2*d3x[2] + 280*sin(x[2])*dx[2]^2*d3x[2]^2 - 280*cos(x[2])*d2x[2]*d3x[2]^2 + 70*cos(x[2])*dx[2]^4*d4x[2] + 420*sin(x[2])*dx[2]^2*d2x[2]*d4x[2] - 210*cos(x[2])*d2x[2]^2*d4x[2] - 280*cos(x[2])*dx[2]*d3x[2]*d4x[2] - 35*sin(x[2])*d4x[2]^2 + 56*sin(x[2])*dx[2]^3*d5x[2] - 168*cos(x[2])*dx[2]*d2x[2]*d5x[2] - 56*sin(x[2])*d3x[2]*d5x[2] - 28*cos(x[2])*dx[2]^2*d6x[2] - 28*sin(x[2])*d2x[2]*d6x[2] - 8*sin(x[2])*dx[2]*d7x[2] + cos(x[2])*d8x[2]) + 9*cos(x[2])*dx[2]*(cos(ξ)*dξ[1]^8 + 28*sin(ξ)*dξ[1]^6*dξ[2] - 210*cos(ξ)*dξ[1]^4*dξ[2]^2 - 420*sin(ξ)*dξ[1]^2*dξ[2]^3 + 105*cos(ξ)*dξ[2]^4 - 56*cos(ξ)*dξ[1]^5*dξ[3] - 560*sin(ξ)*dξ[1]^3*dξ[2]*dξ[3] + 840*cos(ξ)*dξ[1]*dξ[2]^2*dξ[3] + 280*cos(ξ)*dξ[1]^2*dξ[3]^2 + 280*sin(ξ)*dξ[2]*dξ[3]^2 - 70*sin(ξ)*dξ[1]^4*dξ[4] + 420*cos(ξ)*dξ[1]^2*dξ[2]*dξ[4] + 210*sin(ξ)*dξ[2]^2*dξ[4] + 280*sin(ξ)*dξ[1]*dξ[3]*dξ[4] - 35*cos(ξ)*dξ[4]^2 + 56*cos(ξ)*dξ[1]^3*dξ[5] + 168*sin(ξ)*dξ[1]*dξ[2]*dξ[5] - 56*cos(ξ)*dξ[3]*dξ[5] + 28*sin(ξ)*dξ[1]^2*dξ[6] - 28*cos(ξ)*dξ[2]*dξ[6] - 8*cos(ξ)*dξ[1]*dξ[7] - sin(ξ)*dξ[8]) + cos(ξ)*(cos(x[2])*dx[2]^9 + 36*sin(x[2])*dx[2]^7*d2x[2] - 378*cos(x[2])*dx[2]^5*d2x[2]^2 - 1260*sin(x[2])*dx[2]^3*d2x[2]^3 + 945*cos(x[2])*dx[2]*d2x[2]^4 - 84*cos(x[2])*dx[2]^6*d3x[2] - 1260*sin(x[2])*dx[2]^4*d2x[2]*d3x[2] + 3780*cos(x[2])*dx[2]^2*d2x[2]^2*d3x[2] + 1260*sin(x[2])*d2x[2]^3*d3x[2] + 840*cos(x[2])*dx[2]^3*d3x[2]^2 + 2520*sin(x[2])*dx[2]*d2x[2]*d3x[2]^2 - 280*cos(x[2])*d3x[2]^3 - 126*sin(x[2])*dx[2]^5*d4x[2] + 1260*cos(x[2])*dx[2]^3*d2x[2]*d4x[2] + 1890*sin(x[2])*dx[2]*d2x[2]^2*d4x[2] + 1260*sin(x[2])*dx[2]^2*d3x[2]*d4x[2] - 1260*cos(x[2])*d2x[2]*d3x[2]*d4x[2] - 315*cos(x[2])*dx[2]*d4x[2]^2 + 126*cos(x[2])*dx[2]^4*d5x[2] + 756*sin(x[2])*dx[2]^2*d2x[2]*d5x[2] - 378*cos(x[2])*d2x[2]^2*d5x[2] - 504*cos(x[2])*dx[2]*d3x[2]*d5x[2] - 126*sin(x[2])*d4x[2]*d5x[2] + 84*sin(x[2])*dx[2]^3*d6x[2] - 252*cos(x[2])*dx[2]*d2x[2]*d6x[2] - 84*sin(x[2])*d3x[2]*d6x[2] - 36*cos(x[2])*dx[2]^2*d7x[2] - 36*sin(x[2])*d2x[2]*d7x[2] - 9*sin(x[2])*dx[2]*d8x[2] + cos(x[2])*d9x[2]) + sin(x[2])*(-(sin(ξ)*dξ[1]^9) + 36*cos(ξ)*dξ[1]^7*dξ[2] + 378*sin(ξ)*dξ[1]^5*dξ[2]^2 - 1260*cos(ξ)*dξ[1]^3*dξ[2]^3 - 945*sin(ξ)*dξ[1]*dξ[2]^4 + 84*sin(ξ)*dξ[1]^6*dξ[3] - 1260*cos(ξ)*dξ[1]^4*dξ[2]*dξ[3] - 3780*sin(ξ)*dξ[1]^2*dξ[2]^2*dξ[3] + 1260*cos(ξ)*dξ[2]^3*dξ[3] - 840*sin(ξ)*dξ[1]^3*dξ[3]^2 + 2520*cos(ξ)*dξ[1]*dξ[2]*dξ[3]^2 + 280*sin(ξ)*dξ[3]^3 - 126*cos(ξ)*dξ[1]^5*dξ[4] - 1260*sin(ξ)*dξ[1]^3*dξ[2]*dξ[4] + 1890*cos(ξ)*dξ[1]*dξ[2]^2*dξ[4] + 1260*cos(ξ)*dξ[1]^2*dξ[3]*dξ[4] + 1260*sin(ξ)*dξ[2]*dξ[3]*dξ[4] + 315*sin(ξ)*dξ[1]*dξ[4]^2 - 126*sin(ξ)*dξ[1]^4*dξ[5] + 756*cos(ξ)*dξ[1]^2*dξ[2]*dξ[5] + 378*sin(ξ)*dξ[2]^2*dξ[5] + 504*sin(ξ)*dξ[1]*dξ[3]*dξ[5] - 126*cos(ξ)*dξ[4]*dξ[5] + 84*cos(ξ)*dξ[1]^3*dξ[6] + 252*sin(ξ)*dξ[1]*dξ[2]*dξ[6] - 84*cos(ξ)*dξ[3]*dξ[6] + 36*sin(ξ)*dξ[1]^2*dξ[7] - 36*cos(ξ)*dξ[2]*dξ[7] - 9*cos(ξ)*dξ[1]*dξ[8] - sin(ξ)*dξ[9])) + cos(ξ)*sin(x[2])*(280*d3x[1]^3*dρ[3] + 1260*d2x[1]*d3x[1]*dρ[3]*d4x[1] + 315*dx[1]*dρ[3]*d4x[1]^2 + 1260*d2x[1]^3*d3x[1]*dρ[4] + 2520*dx[1]*d2x[1]*d3x[1]^2*dρ[4] + 1890*dx[1]*d2x[1]^2*d4x[1]*dρ[4] + 1260*dx[1]^2*d3x[1]*d4x[1]*dρ[4] + 378*d2x[1]^2*dρ[3]*d5x[1] + 504*dx[1]*d3x[1]*dρ[3]*d5x[1] + 126*dρ[2]*d4x[1]*d5x[1] + 756*dx[1]^2*d2x[1]*dρ[4]*d5x[1] + 945*dx[1]*d2x[1]^4*dρ[5] + 3780*dx[1]^2*d2x[1]^2*d3x[1]*dρ[5] + 840*dx[1]^3*d3x[1]^2*dρ[5] + 1260*dx[1]^3*d2x[1]*d4x[1]*dρ[5] + 126*dx[1]^4*d5x[1]*dρ[5] + 84*dρ[2]*d3x[1]*d6x[1] + 252*dx[1]*d2x[1]*dρ[3]*d6x[1] + 84*dx[1]^3*dρ[4]*d6x[1] + 1260*dx[1]^3*d2x[1]^3*dρ[6] + 1260*dx[1]^4*d2x[1]*d3x[1]*dρ[6] + 126*dx[1]^5*d4x[1]*dρ[6] + 36*d2x[1]*dρ[2]*d7x[1] + 36*dx[1]^2*dρ[3]*d7x[1] + 378*dx[1]^5*d2x[1]^2*dρ[7] + 84*dx[1]^6*d3x[1]*dρ[7] + 9*dx[1]*dρ[2]*d8x[1] + 36*dx[1]^7*d2x[1]*dρ[8] + dρ[1]*d9x[1] + dx[1]^9*dρ[9])

dyH9(x::AbstractVector{Float64}, dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, d4x::AbstractVector{Float64}, d5x::AbstractVector{Float64}, d6x::AbstractVector{Float64}, d7x::AbstractVector{Float64}, d8x::AbstractVector{Float64}, d9x::AbstractVector{Float64}, ξ::Float64, dξ::AbstractVector{Float64}, ρ::Float64, dρ::AbstractVector{Float64}) = 126*(3*d2x[1]^2*dρ[2] + 4*dx[1]*dρ[2]*d3x[1] + 6*dx[1]^2*d2x[1]*dρ[3] + dρ[1]*d4x[1] + dx[1]^4*dρ[4])*(10*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2])*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) + 10*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(-(cos(ξ)*dξ[1]^3) - 3*sin(ξ)*dξ[1]*dξ[2] + cos(ξ)*dξ[3]) + 5*cos(ξ)*dξ[1]*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + 5*cos(x[2])*dx[2]*(sin(ξ)*dξ[1]^4 - 6*cos(ξ)*dξ[1]^2*dξ[2] - 3*sin(ξ)*dξ[2]^2 - 4*sin(ξ)*dξ[1]*dξ[3] + cos(ξ)*dξ[4]) + sin(ξ)*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] - 15*cos(x[2])*dx[2]*d2x[2]^2 - 10*cos(x[2])*dx[2]^2*d3x[2] - 10*sin(x[2])*d2x[2]*d3x[2] - 5*sin(x[2])*dx[2]*d4x[2] + cos(x[2])*d5x[2]) + sin(x[2])*(cos(ξ)*dξ[1]^5 + 10*sin(ξ)*dξ[1]^3*dξ[2] - 15*cos(ξ)*dξ[1]*dξ[2]^2 - 10*cos(ξ)*dξ[1]^2*dξ[3] - 10*sin(ξ)*dξ[2]*dξ[3] - 5*sin(ξ)*dξ[1]*dξ[4] + cos(ξ)*dξ[5])) + 126*(6*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2]) + 4*cos(ξ)*dξ[1]*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) + 4*cos(x[2])*dx[2]*(-(cos(ξ)*dξ[1]^3) - 3*sin(ξ)*dξ[1]*dξ[2] + cos(ξ)*dξ[3]) + sin(ξ)*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + sin(x[2])*(sin(ξ)*dξ[1]^4 - 6*cos(ξ)*dξ[1]^2*dξ[2] - 3*sin(ξ)*dξ[2]^2 - 4*sin(ξ)*dξ[1]*dξ[3] + cos(ξ)*dξ[4]))*(10*d2x[1]*dρ[2]*d3x[1] + 15*dx[1]*d2x[1]^2*dρ[3] + 10*dx[1]^2*d3x[1]*dρ[3] + 5*dx[1]*dρ[2]*d4x[1] + 10*dx[1]^3*d2x[1]*dρ[4] + dρ[1]*d5x[1] + dx[1]^5*dρ[5]) + 84*(3*dx[1]*d2x[1]*dρ[2] + dρ[1]*d3x[1] + dx[1]^3*dρ[3])*(20*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2])*(-(cos(ξ)*dξ[1]^3) - 3*sin(ξ)*dξ[1]*dξ[2] + cos(ξ)*dξ[3]) + 15*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2])*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + 15*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(sin(ξ)*dξ[1]^4 - 6*cos(ξ)*dξ[1]^2*dξ[2] - 3*sin(ξ)*dξ[2]^2 - 4*sin(ξ)*dξ[1]*dξ[3] + cos(ξ)*dξ[4]) + 6*cos(ξ)*dξ[1]*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] - 15*cos(x[2])*dx[2]*d2x[2]^2 - 10*cos(x[2])*dx[2]^2*d3x[2] - 10*sin(x[2])*d2x[2]*d3x[2] - 5*sin(x[2])*dx[2]*d4x[2] + cos(x[2])*d5x[2]) + 6*cos(x[2])*dx[2]*(cos(ξ)*dξ[1]^5 + 10*sin(ξ)*dξ[1]^3*dξ[2] - 15*cos(ξ)*dξ[1]*dξ[2]^2 - 10*cos(ξ)*dξ[1]^2*dξ[3] - 10*sin(ξ)*dξ[2]*dξ[3] - 5*sin(ξ)*dξ[1]*dξ[4] + cos(ξ)*dξ[5]) + sin(ξ)*(-(sin(x[2])*dx[2]^6) + 15*cos(x[2])*dx[2]^4*d2x[2] + 45*sin(x[2])*dx[2]^2*d2x[2]^2 - 15*cos(x[2])*d2x[2]^3 + 20*sin(x[2])*dx[2]^3*d3x[2] - 60*cos(x[2])*dx[2]*d2x[2]*d3x[2] - 10*sin(x[2])*d3x[2]^2 - 15*cos(x[2])*dx[2]^2*d4x[2] - 15*sin(x[2])*d2x[2]*d4x[2] - 6*sin(x[2])*dx[2]*d5x[2] + cos(x[2])*d6x[2]) + sin(x[2])*(-(sin(ξ)*dξ[1]^6) + 15*cos(ξ)*dξ[1]^4*dξ[2] + 45*sin(ξ)*dξ[1]^2*dξ[2]^2 - 15*cos(ξ)*dξ[2]^3 + 20*sin(ξ)*dξ[1]^3*dξ[3] - 60*cos(ξ)*dξ[1]*dξ[2]*dξ[3] - 10*sin(ξ)*dξ[3]^2 - 15*cos(ξ)*dξ[1]^2*dξ[4] - 15*sin(ξ)*dξ[2]*dξ[4] - 6*sin(ξ)*dξ[1]*dξ[5] + cos(ξ)*dξ[6])) + 84*(3*cos(ξ)*dξ[1]*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) + 3*cos(x[2])*dx[2]*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2]) + sin(ξ)*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2]) + sin(x[2])*(-(cos(ξ)*dξ[1]^3) - 3*sin(ξ)*dξ[1]*dξ[2] + cos(ξ)*dξ[3]))*(10*dρ[2]*d3x[1]^2 + 15*d2x[1]^3*dρ[3] + 60*dx[1]*d2x[1]*d3x[1]*dρ[3] + 15*d2x[1]*dρ[2]*d4x[1] + 15*dx[1]^2*dρ[3]*d4x[1] + 45*dx[1]^2*d2x[1]^2*dρ[4] + 20*dx[1]^3*d3x[1]*dρ[4] + 6*dx[1]*dρ[2]*d5x[1] + 15*dx[1]^4*d2x[1]*dρ[5] + dρ[1]*d6x[1] + dx[1]^6*dρ[6]) + 36*(dρ[1]*d2x[1] + dx[1]^2*dρ[2])*(35*(-(cos(ξ)*dξ[1]^3) - 3*sin(ξ)*dξ[1]*dξ[2] + cos(ξ)*dξ[3])*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2]) + 35*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2])*(sin(ξ)*dξ[1]^4 - 6*cos(ξ)*dξ[1]^2*dξ[2] - 3*sin(ξ)*dξ[2]^2 - 4*sin(ξ)*dξ[1]*dξ[3] + cos(ξ)*dξ[4]) + 21*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2])*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] - 15*cos(x[2])*dx[2]*d2x[2]^2 - 10*cos(x[2])*dx[2]^2*d3x[2] - 10*sin(x[2])*d2x[2]*d3x[2] - 5*sin(x[2])*dx[2]*d4x[2] + cos(x[2])*d5x[2]) + 21*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(cos(ξ)*dξ[1]^5 + 10*sin(ξ)*dξ[1]^3*dξ[2] - 15*cos(ξ)*dξ[1]*dξ[2]^2 - 10*cos(ξ)*dξ[1]^2*dξ[3] - 10*sin(ξ)*dξ[2]*dξ[3] - 5*sin(ξ)*dξ[1]*dξ[4] + cos(ξ)*dξ[5]) + 7*cos(ξ)*dξ[1]*(-(sin(x[2])*dx[2]^6) + 15*cos(x[2])*dx[2]^4*d2x[2] + 45*sin(x[2])*dx[2]^2*d2x[2]^2 - 15*cos(x[2])*d2x[2]^3 + 20*sin(x[2])*dx[2]^3*d3x[2] - 60*cos(x[2])*dx[2]*d2x[2]*d3x[2] - 10*sin(x[2])*d3x[2]^2 - 15*cos(x[2])*dx[2]^2*d4x[2] - 15*sin(x[2])*d2x[2]*d4x[2] - 6*sin(x[2])*dx[2]*d5x[2] + cos(x[2])*d6x[2]) + 7*cos(x[2])*dx[2]*(-(sin(ξ)*dξ[1]^6) + 15*cos(ξ)*dξ[1]^4*dξ[2] + 45*sin(ξ)*dξ[1]^2*dξ[2]^2 - 15*cos(ξ)*dξ[2]^3 + 20*sin(ξ)*dξ[1]^3*dξ[3] - 60*cos(ξ)*dξ[1]*dξ[2]*dξ[3] - 10*sin(ξ)*dξ[3]^2 - 15*cos(ξ)*dξ[1]^2*dξ[4] - 15*sin(ξ)*dξ[2]*dξ[4] - 6*sin(ξ)*dξ[1]*dξ[5] + cos(ξ)*dξ[6]) + sin(ξ)*(-(cos(x[2])*dx[2]^7) - 21*sin(x[2])*dx[2]^5*d2x[2] + 105*cos(x[2])*dx[2]^3*d2x[2]^2 + 105*sin(x[2])*dx[2]*d2x[2]^3 + 35*cos(x[2])*dx[2]^4*d3x[2] + 210*sin(x[2])*dx[2]^2*d2x[2]*d3x[2] - 105*cos(x[2])*d2x[2]^2*d3x[2] - 70*cos(x[2])*dx[2]*d3x[2]^2 + 35*sin(x[2])*dx[2]^3*d4x[2] - 105*cos(x[2])*dx[2]*d2x[2]*d4x[2] - 35*sin(x[2])*d3x[2]*d4x[2] - 21*cos(x[2])*dx[2]^2*d5x[2] - 21*sin(x[2])*d2x[2]*d5x[2] - 7*sin(x[2])*dx[2]*d6x[2] + cos(x[2])*d7x[2]) + sin(x[2])*(-(cos(ξ)*dξ[1]^7) - 21*sin(ξ)*dξ[1]^5*dξ[2] + 105*cos(ξ)*dξ[1]^3*dξ[2]^2 + 105*sin(ξ)*dξ[1]*dξ[2]^3 + 35*cos(ξ)*dξ[1]^4*dξ[3] + 210*sin(ξ)*dξ[1]^2*dξ[2]*dξ[3] - 105*cos(ξ)*dξ[2]^2*dξ[3] - 70*cos(ξ)*dξ[1]*dξ[3]^2 + 35*sin(ξ)*dξ[1]^3*dξ[4] - 105*cos(ξ)*dξ[1]*dξ[2]*dξ[4] - 35*sin(ξ)*dξ[3]*dξ[4] - 21*cos(ξ)*dξ[1]^2*dξ[5] - 21*sin(ξ)*dξ[2]*dξ[5] - 7*sin(ξ)*dξ[1]*dξ[6] + cos(ξ)*dξ[7])) + 36*(2*cos(x[2])*cos(ξ)*dx[2]*dξ[1] + sin(ξ)*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2]) + sin(x[2])*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2]))*(105*d2x[1]^2*d3x[1]*dρ[3] + 70*dx[1]*d3x[1]^2*dρ[3] + 35*dρ[2]*d3x[1]*d4x[1] + 105*dx[1]*d2x[1]*dρ[3]*d4x[1] + 105*dx[1]*d2x[1]^3*dρ[4] + 210*dx[1]^2*d2x[1]*d3x[1]*dρ[4] + 35*dx[1]^3*d4x[1]*dρ[4] + 21*d2x[1]*dρ[2]*d5x[1] + 21*dx[1]^2*dρ[3]*d5x[1] + 105*dx[1]^3*d2x[1]^2*dρ[5] + 35*dx[1]^4*d3x[1]*dρ[5] + 7*dx[1]*dρ[2]*d6x[1] + 21*dx[1]^5*d2x[1]*dρ[6] + dρ[1]*d7x[1] + dx[1]^7*dρ[7]) + 9*dx[1]*dρ[1]*(70*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2])*(sin(ξ)*dξ[1]^4 - 6*cos(ξ)*dξ[1]^2*dξ[2] - 3*sin(ξ)*dξ[2]^2 - 4*sin(ξ)*dξ[1]*dξ[3] + cos(ξ)*dξ[4]) + 56*(-(cos(ξ)*dξ[1]^3) - 3*sin(ξ)*dξ[1]*dξ[2] + cos(ξ)*dξ[3])*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] - 15*cos(x[2])*dx[2]*d2x[2]^2 - 10*cos(x[2])*dx[2]^2*d3x[2] - 10*sin(x[2])*d2x[2]*d3x[2] - 5*sin(x[2])*dx[2]*d4x[2] + cos(x[2])*d5x[2]) + 56*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2])*(cos(ξ)*dξ[1]^5 + 10*sin(ξ)*dξ[1]^3*dξ[2] - 15*cos(ξ)*dξ[1]*dξ[2]^2 - 10*cos(ξ)*dξ[1]^2*dξ[3] - 10*sin(ξ)*dξ[2]*dξ[3] - 5*sin(ξ)*dξ[1]*dξ[4] + cos(ξ)*dξ[5]) + 28*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2])*(-(sin(x[2])*dx[2]^6) + 15*cos(x[2])*dx[2]^4*d2x[2] + 45*sin(x[2])*dx[2]^2*d2x[2]^2 - 15*cos(x[2])*d2x[2]^3 + 20*sin(x[2])*dx[2]^3*d3x[2] - 60*cos(x[2])*dx[2]*d2x[2]*d3x[2] - 10*sin(x[2])*d3x[2]^2 - 15*cos(x[2])*dx[2]^2*d4x[2] - 15*sin(x[2])*d2x[2]*d4x[2] - 6*sin(x[2])*dx[2]*d5x[2] + cos(x[2])*d6x[2]) + 28*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(-(sin(ξ)*dξ[1]^6) + 15*cos(ξ)*dξ[1]^4*dξ[2] + 45*sin(ξ)*dξ[1]^2*dξ[2]^2 - 15*cos(ξ)*dξ[2]^3 + 20*sin(ξ)*dξ[1]^3*dξ[3] - 60*cos(ξ)*dξ[1]*dξ[2]*dξ[3] - 10*sin(ξ)*dξ[3]^2 - 15*cos(ξ)*dξ[1]^2*dξ[4] - 15*sin(ξ)*dξ[2]*dξ[4] - 6*sin(ξ)*dξ[1]*dξ[5] + cos(ξ)*dξ[6]) + 8*cos(ξ)*dξ[1]*(-(cos(x[2])*dx[2]^7) - 21*sin(x[2])*dx[2]^5*d2x[2] + 105*cos(x[2])*dx[2]^3*d2x[2]^2 + 105*sin(x[2])*dx[2]*d2x[2]^3 + 35*cos(x[2])*dx[2]^4*d3x[2] + 210*sin(x[2])*dx[2]^2*d2x[2]*d3x[2] - 105*cos(x[2])*d2x[2]^2*d3x[2] - 70*cos(x[2])*dx[2]*d3x[2]^2 + 35*sin(x[2])*dx[2]^3*d4x[2] - 105*cos(x[2])*dx[2]*d2x[2]*d4x[2] - 35*sin(x[2])*d3x[2]*d4x[2] - 21*cos(x[2])*dx[2]^2*d5x[2] - 21*sin(x[2])*d2x[2]*d5x[2] - 7*sin(x[2])*dx[2]*d6x[2] + cos(x[2])*d7x[2]) + 8*cos(x[2])*dx[2]*(-(cos(ξ)*dξ[1]^7) - 21*sin(ξ)*dξ[1]^5*dξ[2] + 105*cos(ξ)*dξ[1]^3*dξ[2]^2 + 105*sin(ξ)*dξ[1]*dξ[2]^3 + 35*cos(ξ)*dξ[1]^4*dξ[3] + 210*sin(ξ)*dξ[1]^2*dξ[2]*dξ[3] - 105*cos(ξ)*dξ[2]^2*dξ[3] - 70*cos(ξ)*dξ[1]*dξ[3]^2 + 35*sin(ξ)*dξ[1]^3*dξ[4] - 105*cos(ξ)*dξ[1]*dξ[2]*dξ[4] - 35*sin(ξ)*dξ[3]*dξ[4] - 21*cos(ξ)*dξ[1]^2*dξ[5] - 21*sin(ξ)*dξ[2]*dξ[5] - 7*sin(ξ)*dξ[1]*dξ[6] + cos(ξ)*dξ[7]) + sin(ξ)*(sin(x[2])*dx[2]^8 - 28*cos(x[2])*dx[2]^6*d2x[2] - 210*sin(x[2])*dx[2]^4*d2x[2]^2 + 420*cos(x[2])*dx[2]^2*d2x[2]^3 + 105*sin(x[2])*d2x[2]^4 - 56*sin(x[2])*dx[2]^5*d3x[2] + 560*cos(x[2])*dx[2]^3*d2x[2]*d3x[2] + 840*sin(x[2])*dx[2]*d2x[2]^2*d3x[2] + 280*sin(x[2])*dx[2]^2*d3x[2]^2 - 280*cos(x[2])*d2x[2]*d3x[2]^2 + 70*cos(x[2])*dx[2]^4*d4x[2] + 420*sin(x[2])*dx[2]^2*d2x[2]*d4x[2] - 210*cos(x[2])*d2x[2]^2*d4x[2] - 280*cos(x[2])*dx[2]*d3x[2]*d4x[2] - 35*sin(x[2])*d4x[2]^2 + 56*sin(x[2])*dx[2]^3*d5x[2] - 168*cos(x[2])*dx[2]*d2x[2]*d5x[2] - 56*sin(x[2])*d3x[2]*d5x[2] - 28*cos(x[2])*dx[2]^2*d6x[2] - 28*sin(x[2])*d2x[2]*d6x[2] - 8*sin(x[2])*dx[2]*d7x[2] + cos(x[2])*d8x[2]) + sin(x[2])*(sin(ξ)*dξ[1]^8 - 28*cos(ξ)*dξ[1]^6*dξ[2] - 210*sin(ξ)*dξ[1]^4*dξ[2]^2 + 420*cos(ξ)*dξ[1]^2*dξ[2]^3 + 105*sin(ξ)*dξ[2]^4 - 56*sin(ξ)*dξ[1]^5*dξ[3] + 560*cos(ξ)*dξ[1]^3*dξ[2]*dξ[3] + 840*sin(ξ)*dξ[1]*dξ[2]^2*dξ[3] + 280*sin(ξ)*dξ[1]^2*dξ[3]^2 - 280*cos(ξ)*dξ[2]*dξ[3]^2 + 70*cos(ξ)*dξ[1]^4*dξ[4] + 420*sin(ξ)*dξ[1]^2*dξ[2]*dξ[4] - 210*cos(ξ)*dξ[2]^2*dξ[4] - 280*cos(ξ)*dξ[1]*dξ[3]*dξ[4] - 35*sin(ξ)*dξ[4]^2 + 56*sin(ξ)*dξ[1]^3*dξ[5] - 168*cos(ξ)*dξ[1]*dξ[2]*dξ[5] - 56*sin(ξ)*dξ[3]*dξ[5] - 28*cos(ξ)*dξ[1]^2*dξ[6] - 28*sin(ξ)*dξ[2]*dξ[6] - 8*sin(ξ)*dξ[1]*dξ[7] + cos(ξ)*dξ[8])) + 9*(cos(x[2])*sin(ξ)*dx[2] + cos(ξ)*sin(x[2])*dξ[1])*(280*d2x[1]*d3x[1]^2*dρ[3] + 210*d2x[1]^2*dρ[3]*d4x[1] + 280*dx[1]*d3x[1]*dρ[3]*d4x[1] + 35*dρ[2]*d4x[1]^2 + 105*d2x[1]^4*dρ[4] + 840*dx[1]*d2x[1]^2*d3x[1]*dρ[4] + 280*dx[1]^2*d3x[1]^2*dρ[4] + 420*dx[1]^2*d2x[1]*d4x[1]*dρ[4] + 56*dρ[2]*d3x[1]*d5x[1] + 168*dx[1]*d2x[1]*dρ[3]*d5x[1] + 56*dx[1]^3*dρ[4]*d5x[1] + 420*dx[1]^2*d2x[1]^3*dρ[5] + 560*dx[1]^3*d2x[1]*d3x[1]*dρ[5] + 70*dx[1]^4*d4x[1]*dρ[5] + 28*d2x[1]*dρ[2]*d6x[1] + 28*dx[1]^2*dρ[3]*d6x[1] + 210*dx[1]^4*d2x[1]^2*dρ[6] + 56*dx[1]^5*d3x[1]*dρ[6] + 8*dx[1]*dρ[2]*d7x[1] + 28*dx[1]^6*d2x[1]*dρ[7] + dρ[1]*d8x[1] + dx[1]^8*dρ[8]) + ρ*(126*(sin(ξ)*dξ[1]^4 - 6*cos(ξ)*dξ[1]^2*dξ[2] - 3*sin(ξ)*dξ[2]^2 - 4*sin(ξ)*dξ[1]*dξ[3] + cos(ξ)*dξ[4])*(cos(x[2])*dx[2]^5 + 10*sin(x[2])*dx[2]^3*d2x[2] - 15*cos(x[2])*dx[2]*d2x[2]^2 - 10*cos(x[2])*dx[2]^2*d3x[2] - 10*sin(x[2])*d2x[2]*d3x[2] - 5*sin(x[2])*dx[2]*d4x[2] + cos(x[2])*d5x[2]) + 126*(sin(x[2])*dx[2]^4 - 6*cos(x[2])*dx[2]^2*d2x[2] - 3*sin(x[2])*d2x[2]^2 - 4*sin(x[2])*dx[2]*d3x[2] + cos(x[2])*d4x[2])*(cos(ξ)*dξ[1]^5 + 10*sin(ξ)*dξ[1]^3*dξ[2] - 15*cos(ξ)*dξ[1]*dξ[2]^2 - 10*cos(ξ)*dξ[1]^2*dξ[3] - 10*sin(ξ)*dξ[2]*dξ[3] - 5*sin(ξ)*dξ[1]*dξ[4] + cos(ξ)*dξ[5]) + 84*(-(cos(ξ)*dξ[1]^3) - 3*sin(ξ)*dξ[1]*dξ[2] + cos(ξ)*dξ[3])*(-(sin(x[2])*dx[2]^6) + 15*cos(x[2])*dx[2]^4*d2x[2] + 45*sin(x[2])*dx[2]^2*d2x[2]^2 - 15*cos(x[2])*d2x[2]^3 + 20*sin(x[2])*dx[2]^3*d3x[2] - 60*cos(x[2])*dx[2]*d2x[2]*d3x[2] - 10*sin(x[2])*d3x[2]^2 - 15*cos(x[2])*dx[2]^2*d4x[2] - 15*sin(x[2])*d2x[2]*d4x[2] - 6*sin(x[2])*dx[2]*d5x[2] + cos(x[2])*d6x[2]) + 84*(-(cos(x[2])*dx[2]^3) - 3*sin(x[2])*dx[2]*d2x[2] + cos(x[2])*d3x[2])*(-(sin(ξ)*dξ[1]^6) + 15*cos(ξ)*dξ[1]^4*dξ[2] + 45*sin(ξ)*dξ[1]^2*dξ[2]^2 - 15*cos(ξ)*dξ[2]^3 + 20*sin(ξ)*dξ[1]^3*dξ[3] - 60*cos(ξ)*dξ[1]*dξ[2]*dξ[3] - 10*sin(ξ)*dξ[3]^2 - 15*cos(ξ)*dξ[1]^2*dξ[4] - 15*sin(ξ)*dξ[2]*dξ[4] - 6*sin(ξ)*dξ[1]*dξ[5] + cos(ξ)*dξ[6]) + 36*(-(sin(ξ)*dξ[1]^2) + cos(ξ)*dξ[2])*(-(cos(x[2])*dx[2]^7) - 21*sin(x[2])*dx[2]^5*d2x[2] + 105*cos(x[2])*dx[2]^3*d2x[2]^2 + 105*sin(x[2])*dx[2]*d2x[2]^3 + 35*cos(x[2])*dx[2]^4*d3x[2] + 210*sin(x[2])*dx[2]^2*d2x[2]*d3x[2] - 105*cos(x[2])*d2x[2]^2*d3x[2] - 70*cos(x[2])*dx[2]*d3x[2]^2 + 35*sin(x[2])*dx[2]^3*d4x[2] - 105*cos(x[2])*dx[2]*d2x[2]*d4x[2] - 35*sin(x[2])*d3x[2]*d4x[2] - 21*cos(x[2])*dx[2]^2*d5x[2] - 21*sin(x[2])*d2x[2]*d5x[2] - 7*sin(x[2])*dx[2]*d6x[2] + cos(x[2])*d7x[2]) + 36*(-(sin(x[2])*dx[2]^2) + cos(x[2])*d2x[2])*(-(cos(ξ)*dξ[1]^7) - 21*sin(ξ)*dξ[1]^5*dξ[2] + 105*cos(ξ)*dξ[1]^3*dξ[2]^2 + 105*sin(ξ)*dξ[1]*dξ[2]^3 + 35*cos(ξ)*dξ[1]^4*dξ[3] + 210*sin(ξ)*dξ[1]^2*dξ[2]*dξ[3] - 105*cos(ξ)*dξ[2]^2*dξ[3] - 70*cos(ξ)*dξ[1]*dξ[3]^2 + 35*sin(ξ)*dξ[1]^3*dξ[4] - 105*cos(ξ)*dξ[1]*dξ[2]*dξ[4] - 35*sin(ξ)*dξ[3]*dξ[4] - 21*cos(ξ)*dξ[1]^2*dξ[5] - 21*sin(ξ)*dξ[2]*dξ[5] - 7*sin(ξ)*dξ[1]*dξ[6] + cos(ξ)*dξ[7]) + 9*cos(ξ)*dξ[1]*(sin(x[2])*dx[2]^8 - 28*cos(x[2])*dx[2]^6*d2x[2] - 210*sin(x[2])*dx[2]^4*d2x[2]^2 + 420*cos(x[2])*dx[2]^2*d2x[2]^3 + 105*sin(x[2])*d2x[2]^4 - 56*sin(x[2])*dx[2]^5*d3x[2] + 560*cos(x[2])*dx[2]^3*d2x[2]*d3x[2] + 840*sin(x[2])*dx[2]*d2x[2]^2*d3x[2] + 280*sin(x[2])*dx[2]^2*d3x[2]^2 - 280*cos(x[2])*d2x[2]*d3x[2]^2 + 70*cos(x[2])*dx[2]^4*d4x[2] + 420*sin(x[2])*dx[2]^2*d2x[2]*d4x[2] - 210*cos(x[2])*d2x[2]^2*d4x[2] - 280*cos(x[2])*dx[2]*d3x[2]*d4x[2] - 35*sin(x[2])*d4x[2]^2 + 56*sin(x[2])*dx[2]^3*d5x[2] - 168*cos(x[2])*dx[2]*d2x[2]*d5x[2] - 56*sin(x[2])*d3x[2]*d5x[2] - 28*cos(x[2])*dx[2]^2*d6x[2] - 28*sin(x[2])*d2x[2]*d6x[2] - 8*sin(x[2])*dx[2]*d7x[2] + cos(x[2])*d8x[2]) + 9*cos(x[2])*dx[2]*(sin(ξ)*dξ[1]^8 - 28*cos(ξ)*dξ[1]^6*dξ[2] - 210*sin(ξ)*dξ[1]^4*dξ[2]^2 + 420*cos(ξ)*dξ[1]^2*dξ[2]^3 + 105*sin(ξ)*dξ[2]^4 - 56*sin(ξ)*dξ[1]^5*dξ[3] + 560*cos(ξ)*dξ[1]^3*dξ[2]*dξ[3] + 840*sin(ξ)*dξ[1]*dξ[2]^2*dξ[3] + 280*sin(ξ)*dξ[1]^2*dξ[3]^2 - 280*cos(ξ)*dξ[2]*dξ[3]^2 + 70*cos(ξ)*dξ[1]^4*dξ[4] + 420*sin(ξ)*dξ[1]^2*dξ[2]*dξ[4] - 210*cos(ξ)*dξ[2]^2*dξ[4] - 280*cos(ξ)*dξ[1]*dξ[3]*dξ[4] - 35*sin(ξ)*dξ[4]^2 + 56*sin(ξ)*dξ[1]^3*dξ[5] - 168*cos(ξ)*dξ[1]*dξ[2]*dξ[5] - 56*sin(ξ)*dξ[3]*dξ[5] - 28*cos(ξ)*dξ[1]^2*dξ[6] - 28*sin(ξ)*dξ[2]*dξ[6] - 8*sin(ξ)*dξ[1]*dξ[7] + cos(ξ)*dξ[8]) + sin(ξ)*(cos(x[2])*dx[2]^9 + 36*sin(x[2])*dx[2]^7*d2x[2] - 378*cos(x[2])*dx[2]^5*d2x[2]^2 - 1260*sin(x[2])*dx[2]^3*d2x[2]^3 + 945*cos(x[2])*dx[2]*d2x[2]^4 - 84*cos(x[2])*dx[2]^6*d3x[2] - 1260*sin(x[2])*dx[2]^4*d2x[2]*d3x[2] + 3780*cos(x[2])*dx[2]^2*d2x[2]^2*d3x[2] + 1260*sin(x[2])*d2x[2]^3*d3x[2] + 840*cos(x[2])*dx[2]^3*d3x[2]^2 + 2520*sin(x[2])*dx[2]*d2x[2]*d3x[2]^2 - 280*cos(x[2])*d3x[2]^3 - 126*sin(x[2])*dx[2]^5*d4x[2] + 1260*cos(x[2])*dx[2]^3*d2x[2]*d4x[2] + 1890*sin(x[2])*dx[2]*d2x[2]^2*d4x[2] + 1260*sin(x[2])*dx[2]^2*d3x[2]*d4x[2] - 1260*cos(x[2])*d2x[2]*d3x[2]*d4x[2] - 315*cos(x[2])*dx[2]*d4x[2]^2 + 126*cos(x[2])*dx[2]^4*d5x[2] + 756*sin(x[2])*dx[2]^2*d2x[2]*d5x[2] - 378*cos(x[2])*d2x[2]^2*d5x[2] - 504*cos(x[2])*dx[2]*d3x[2]*d5x[2] - 126*sin(x[2])*d4x[2]*d5x[2] + 84*sin(x[2])*dx[2]^3*d6x[2] - 252*cos(x[2])*dx[2]*d2x[2]*d6x[2] - 84*sin(x[2])*d3x[2]*d6x[2] - 36*cos(x[2])*dx[2]^2*d7x[2] - 36*sin(x[2])*d2x[2]*d7x[2] - 9*sin(x[2])*dx[2]*d8x[2] + cos(x[2])*d9x[2]) + sin(x[2])*(cos(ξ)*dξ[1]^9 + 36*sin(ξ)*dξ[1]^7*dξ[2] - 378*cos(ξ)*dξ[1]^5*dξ[2]^2 - 1260*sin(ξ)*dξ[1]^3*dξ[2]^3 + 945*cos(ξ)*dξ[1]*dξ[2]^4 - 84*cos(ξ)*dξ[1]^6*dξ[3] - 1260*sin(ξ)*dξ[1]^4*dξ[2]*dξ[3] + 3780*cos(ξ)*dξ[1]^2*dξ[2]^2*dξ[3] + 1260*sin(ξ)*dξ[2]^3*dξ[3] + 840*cos(ξ)*dξ[1]^3*dξ[3]^2 + 2520*sin(ξ)*dξ[1]*dξ[2]*dξ[3]^2 - 280*cos(ξ)*dξ[3]^3 - 126*sin(ξ)*dξ[1]^5*dξ[4] + 1260*cos(ξ)*dξ[1]^3*dξ[2]*dξ[4] + 1890*sin(ξ)*dξ[1]*dξ[2]^2*dξ[4] + 1260*sin(ξ)*dξ[1]^2*dξ[3]*dξ[4] - 1260*cos(ξ)*dξ[2]*dξ[3]*dξ[4] - 315*cos(ξ)*dξ[1]*dξ[4]^2 + 126*cos(ξ)*dξ[1]^4*dξ[5] + 756*sin(ξ)*dξ[1]^2*dξ[2]*dξ[5] - 378*cos(ξ)*dξ[2]^2*dξ[5] - 504*cos(ξ)*dξ[1]*dξ[3]*dξ[5] - 126*sin(ξ)*dξ[4]*dξ[5] + 84*sin(ξ)*dξ[1]^3*dξ[6] - 252*cos(ξ)*dξ[1]*dξ[2]*dξ[6] - 84*sin(ξ)*dξ[3]*dξ[6] - 36*cos(ξ)*dξ[1]^2*dξ[7] - 36*sin(ξ)*dξ[2]*dξ[7] - 9*sin(ξ)*dξ[1]*dξ[8] + cos(ξ)*dξ[9])) + sin(x[2])*sin(ξ)*(280*d3x[1]^3*dρ[3] + 1260*d2x[1]*d3x[1]*dρ[3]*d4x[1] + 315*dx[1]*dρ[3]*d4x[1]^2 + 1260*d2x[1]^3*d3x[1]*dρ[4] + 2520*dx[1]*d2x[1]*d3x[1]^2*dρ[4] + 1890*dx[1]*d2x[1]^2*d4x[1]*dρ[4] + 1260*dx[1]^2*d3x[1]*d4x[1]*dρ[4] + 378*d2x[1]^2*dρ[3]*d5x[1] + 504*dx[1]*d3x[1]*dρ[3]*d5x[1] + 126*dρ[2]*d4x[1]*d5x[1] + 756*dx[1]^2*d2x[1]*dρ[4]*d5x[1] + 945*dx[1]*d2x[1]^4*dρ[5] + 3780*dx[1]^2*d2x[1]^2*d3x[1]*dρ[5] + 840*dx[1]^3*d3x[1]^2*dρ[5] + 1260*dx[1]^3*d2x[1]*d4x[1]*dρ[5] + 126*dx[1]^4*d5x[1]*dρ[5] + 84*dρ[2]*d3x[1]*d6x[1] + 252*dx[1]*d2x[1]*dρ[3]*d6x[1] + 84*dx[1]^3*dρ[4]*d6x[1] + 1260*dx[1]^3*d2x[1]^3*dρ[6] + 1260*dx[1]^4*d2x[1]*d3x[1]*dρ[6] + 126*dx[1]^5*d4x[1]*dρ[6] + 36*d2x[1]*dρ[2]*d7x[1] + 36*dx[1]^2*dρ[3]*d7x[1] + 378*dx[1]^5*d2x[1]^2*dρ[7] + 84*dx[1]^6*d3x[1]*dρ[7] + 9*dx[1]*dρ[2]*d8x[1] + 36*dx[1]^7*d2x[1]*dρ[8] + dρ[1]*d9x[1] + dx[1]^9*dρ[9])

dzH9(x::AbstractVector{Float64}, dx::AbstractVector{Float64}, d2x::AbstractVector{Float64}, d3x::AbstractVector{Float64}, d4x::AbstractVector{Float64}, d5x::AbstractVector{Float64}, d6x::AbstractVector{Float64}, d7x::AbstractVector{Float64}, d8x::AbstractVector{Float64}, d9x::AbstractVector{Float64}) = 126*(cos(x[2])*dx[2]^4 + 6*sin(x[2])*dx[2]^2*d2x[2] - 3*cos(x[2])*d2x[2]^2 - 4*cos(x[2])*dx[2]*d3x[2] - sin(x[2])*d4x[2])*d5x[1] + 126*d4x[1]*(-(sin(x[2])*dx[2]^5) + 10*cos(x[2])*dx[2]^3*d2x[2] + 15*sin(x[2])*dx[2]*d2x[2]^2 + 10*sin(x[2])*dx[2]^2*d3x[2] - 10*cos(x[2])*d2x[2]*d3x[2] - 5*cos(x[2])*dx[2]*d4x[2] - sin(x[2])*d5x[2]) + 84*(sin(x[2])*dx[2]^3 - 3*cos(x[2])*dx[2]*d2x[2] - sin(x[2])*d3x[2])*d6x[1] + 84*d3x[1]*(-(cos(x[2])*dx[2]^6) - 15*sin(x[2])*dx[2]^4*d2x[2] + 45*cos(x[2])*dx[2]^2*d2x[2]^2 + 15*sin(x[2])*d2x[2]^3 + 20*cos(x[2])*dx[2]^3*d3x[2] + 60*sin(x[2])*dx[2]*d2x[2]*d3x[2] - 10*cos(x[2])*d3x[2]^2 + 15*sin(x[2])*dx[2]^2*d4x[2] - 15*cos(x[2])*d2x[2]*d4x[2] - 6*cos(x[2])*dx[2]*d5x[2] - sin(x[2])*d6x[2]) + 36*(-(cos(x[2])*dx[2]^2) - sin(x[2])*d2x[2])*d7x[1] + 36*d2x[1]*(sin(x[2])*dx[2]^7 - 21*cos(x[2])*dx[2]^5*d2x[2] - 105*sin(x[2])*dx[2]^3*d2x[2]^2 + 105*cos(x[2])*dx[2]*d2x[2]^3 - 35*sin(x[2])*dx[2]^4*d3x[2] + 210*cos(x[2])*dx[2]^2*d2x[2]*d3x[2] + 105*sin(x[2])*d2x[2]^2*d3x[2] + 70*sin(x[2])*dx[2]*d3x[2]^2 + 35*cos(x[2])*dx[2]^3*d4x[2] + 105*sin(x[2])*dx[2]*d2x[2]*d4x[2] - 35*cos(x[2])*d3x[2]*d4x[2] + 21*sin(x[2])*dx[2]^2*d5x[2] - 21*cos(x[2])*d2x[2]*d5x[2] - 7*cos(x[2])*dx[2]*d6x[2] - sin(x[2])*d7x[2]) - 9*sin(x[2])*dx[2]*d8x[1] + 9*dx[1]*(cos(x[2])*dx[2]^8 + 28*sin(x[2])*dx[2]^6*d2x[2] - 210*cos(x[2])*dx[2]^4*d2x[2]^2 - 420*sin(x[2])*dx[2]^2*d2x[2]^3 + 105*cos(x[2])*d2x[2]^4 - 56*cos(x[2])*dx[2]^5*d3x[2] - 560*sin(x[2])*dx[2]^3*d2x[2]*d3x[2] + 840*cos(x[2])*dx[2]*d2x[2]^2*d3x[2] + 280*cos(x[2])*dx[2]^2*d3x[2]^2 + 280*sin(x[2])*d2x[2]*d3x[2]^2 - 70*sin(x[2])*dx[2]^4*d4x[2] + 420*cos(x[2])*dx[2]^2*d2x[2]*d4x[2] + 210*sin(x[2])*d2x[2]^2*d4x[2] + 280*sin(x[2])*dx[2]*d3x[2]*d4x[2] - 35*cos(x[2])*d4x[2]^2 + 56*cos(x[2])*dx[2]^3*d5x[2] + 168*sin(x[2])*dx[2]*d2x[2]*d5x[2] - 56*cos(x[2])*d3x[2]*d5x[2] + 28*sin(x[2])*dx[2]^2*d6x[2] - 28*cos(x[2])*d2x[2]*d6x[2] - 8*cos(x[2])*dx[2]*d7x[2] - sin(x[2])*d8x[2]) + cos(x[2])*d9x[1] + (-M + x[1])*(-(sin(x[2])*dx[2]^9) + 36*cos(x[2])*dx[2]^7*d2x[2] + 378*sin(x[2])*dx[2]^5*d2x[2]^2 - 1260*cos(x[2])*dx[2]^3*d2x[2]^3 - 945*sin(x[2])*dx[2]*d2x[2]^4 + 84*sin(x[2])*dx[2]^6*d3x[2] - 1260*cos(x[2])*dx[2]^4*d2x[2]*d3x[2] - 3780*sin(x[2])*dx[2]^2*d2x[2]^2*d3x[2] + 1260*cos(x[2])*d2x[2]^3*d3x[2] - 840*sin(x[2])*dx[2]^3*d3x[2]^2 + 2520*cos(x[2])*dx[2]*d2x[2]*d3x[2]^2 + 280*sin(x[2])*d3x[2]^3 - 126*cos(x[2])*dx[2]^5*d4x[2] - 1260*sin(x[2])*dx[2]^3*d2x[2]*d4x[2] + 1890*cos(x[2])*dx[2]*d2x[2]^2*d4x[2] + 1260*cos(x[2])*dx[2]^2*d3x[2]*d4x[2] + 1260*sin(x[2])*d2x[2]*d3x[2]*d4x[2] + 315*sin(x[2])*dx[2]*d4x[2]^2 - 126*sin(x[2])*dx[2]^4*d5x[2] + 756*cos(x[2])*dx[2]^2*d2x[2]*d5x[2] + 378*sin(x[2])*d2x[2]^2*d5x[2] + 504*sin(x[2])*dx[2]*d3x[2]*d5x[2] - 126*cos(x[2])*d4x[2]*d5x[2] + 84*cos(x[2])*dx[2]^3*d6x[2] + 252*sin(x[2])*dx[2]*d2x[2]*d6x[2] - 84*cos(x[2])*d3x[2]*d6x[2] + 36*sin(x[2])*dx[2]^2*d7x[2] - 36*cos(x[2])*d2x[2]*d7x[2] - 9*cos(x[2])*dx[2]*d8x[2] - sin(x[2])*d9x[2])



function compute_harmonic_derivs!(xBL::AbstractVector{Float64}, dxBL::AbstractVector{Float64}, d2xBL::AbstractVector{Float64}, d3xBL::AbstractVector{Float64}, d4xBL::AbstractVector{Float64}, d5xBL::AbstractVector{Float64}, d6xBL::AbstractVector{Float64}, d7xBL::AbstractVector{Float64}, d8xBL::AbstractVector{Float64}, d9xBL::AbstractVector{Float64}, xH::AbstractVector{Float64}, dxH::AbstractVector{Float64}, d2xH::AbstractVector{Float64}, d3xH::AbstractVector{Float64}, d4xH::AbstractVector{Float64}, d5xH::AbstractVector{Float64}, d6xH::AbstractVector{Float64}, d7xH::AbstractVector{Float64}, d8xH::AbstractVector{Float64}, d9xH::AbstractVector{Float64}, a::Float64)

    # inner and outer horizon
    rminus = M - sqrt(M^2 - a^2)
    rplus = M + sqrt(M^2 - a^2)
    r = xBL[1]

    # compute derivatives of ρ
    ρ = HarmonicCoordDerivs.ρ(r, a)
    dρ = zeros(9)
    compute_ρ_derivs!(dρ, r, a)

    # compute derivatives of ξ
    Φ = HarmonicCoordDerivs.Φ(r, a, rminus, rplus)
    dΦ = zeros(9)
    compute_Φ_derivs!(dΦ, r, a, rminus, rplus)

    ξ = HarmonicCoordDerivs.ξ(xBL, Φ)
    dξ = zeros(9)
    compute_ξ_derivs!(dξ, dxBL, d2xBL, d3xBL, d4xBL, d5xBL, d6xBL, d7xBL, d8xBL, d9xBL, dΦ)

    # compute harmonic coordinate derivatives
    xH[1] = HarmonicCoordDerivs.xH(xBL, ξ, ρ)
    xH[2] = HarmonicCoordDerivs.yH(xBL, ξ, ρ)
    xH[3] = HarmonicCoordDerivs.zH(xBL, r)

    dxH[1] = dxH1(xBL, dxBL, ξ, dξ, ρ, dρ)
    dxH[2] = dyH1(xBL, dxBL, ξ, dξ, ρ, dρ)
    dxH[3] = dzH1(xBL, dxBL)

    d2xH[1] = dxH2(xBL, dxBL, d2xBL, ξ, dξ, ρ, dρ)
    d2xH[2] = dyH2(xBL, dxBL, d2xBL, ξ, dξ, ρ, dρ)
    d2xH[3] = dzH2(xBL, dxBL, d2xBL)

    d3xH[1] = dxH3(xBL, dxBL, d2xBL, d3xBL, ξ, dξ, ρ, dρ)
    d3xH[2] = dyH3(xBL, dxBL, d2xBL, d3xBL, ξ, dξ, ρ, dρ)
    d3xH[3] = dzH3(xBL, dxBL, d2xBL, d3xBL)

    d4xH[1] = dxH4(xBL, dxBL, d2xBL, d3xBL, d4xBL, ξ, dξ, ρ, dρ)
    d4xH[2] = dyH4(xBL, dxBL, d2xBL, d3xBL, d4xBL, ξ, dξ, ρ, dρ)
    d4xH[3] = dzH4(xBL, dxBL, d2xBL, d3xBL, d4xBL)

    d5xH[1] = dxH5(xBL, dxBL, d2xBL, d3xBL, d4xBL, d5xBL, ξ, dξ, ρ, dρ)
    d5xH[2] = dyH5(xBL, dxBL, d2xBL, d3xBL, d4xBL, d5xBL, ξ, dξ, ρ, dρ)
    d5xH[3] = dzH5(xBL, dxBL, d2xBL, d3xBL, d4xBL, d5xBL)

    d6xH[1] = dxH6(xBL, dxBL, d2xBL, d3xBL, d4xBL, d5xBL, d6xBL, ξ, dξ, ρ, dρ)
    d6xH[2] = dyH6(xBL, dxBL, d2xBL, d3xBL, d4xBL, d5xBL, d6xBL, ξ, dξ, ρ, dρ)
    d6xH[3] = dzH6(xBL, dxBL, d2xBL, d3xBL, d4xBL, d5xBL, d6xBL)

    d7xH[1] = dxH7(xBL, dxBL, d2xBL, d3xBL, d4xBL, d5xBL, d6xBL, d7xBL, ξ, dξ, ρ, dρ)
    d7xH[2] = dyH7(xBL, dxBL, d2xBL, d3xBL, d4xBL, d5xBL, d6xBL, d7xBL, ξ, dξ, ρ, dρ)
    d7xH[3] = dzH7(xBL, dxBL, d2xBL, d3xBL, d4xBL, d5xBL, d6xBL, d7xBL)

    d8xH[1] = dxH8(xBL, dxBL, d2xBL, d3xBL, d4xBL, d5xBL, d6xBL, d7xBL, d8xBL, ξ, dξ, ρ, dρ)
    d8xH[2] = dyH8(xBL, dxBL, d2xBL, d3xBL, d4xBL, d5xBL, d6xBL, d7xBL, d8xBL, ξ, dξ, ρ, dρ)
    d8xH[3] = dzH8(xBL, dxBL, d2xBL, d3xBL, d4xBL, d5xBL, d6xBL, d7xBL, d8xBL)

    d9xH[1] = dxH9(xBL, dxBL, d2xBL, d3xBL, d4xBL, d5xBL, d6xBL, d7xBL, d8xBL, d9xBL, ξ, dξ, ρ, dρ)
    d9xH[2] = dyH9(xBL, dxBL, d2xBL, d3xBL, d4xBL, d5xBL, d6xBL, d7xBL, d8xBL, d9xBL, ξ, dξ, ρ, dρ)
    d9xH[3] = dzH9(xBL, dxBL, d2xBL, d3xBL, d4xBL, d5xBL, d6xBL, d7xBL, d8xBL, d9xBL)
end

end

# rD(n) = 3 + n/2
# thetaD(n) = π/3 + n/3
# phiD(n) = π/8 + n/4

# xBL = [rD(0), thetaD(0), phiD(0)];
# dxBL = [rD(1), thetaD(1), phiD(1)];
# d2xBL = [rD(2), thetaD(2), phiD(2)];
# d3xBL = [rD(3), thetaD(3), phiD(3)];
# d4xBL = [rD(4), thetaD(4), phiD(4)];
# d5xBL = [rD(5), thetaD(5), phiD(5)];
# d6xBL = [rD(6), thetaD(6), phiD(6)];
# d7xBL = [rD(7), thetaD(7), phiD(7)];
# d8xBL = [rD(8), thetaD(8), phiD(8)];
# d9xBL = [rD(9), thetaD(9), phiD(9)];

# xH = zeros(3);
# dxH = zeros(3);
# d2xH = zeros(3);
# d3xH = zeros(3);
# d4xH = zeros(3);
# d5xH = zeros(3);
# d6xH = zeros(3);
# d7xH = zeros(3);
# d8xH = zeros(3);
# d9xH = zeros(3);
# a = 0.5;

# xHDerivs_MMA = [1.6335363106328038,3.620782580072471,2.659017384140555,-25.47758391755284,-198.2293766241397,-604.3710144278604,-2103.6852042846485,50196.72502762086,-302035.78734692995,1.3420726353825215*1e7]
# yHDerivs_MMA = [0.7204575781017,2.603785045887765,9.433240825620514,14.424221724291959,-25.315039292192928,-1100.1181048031665,-2817.6511045174466,-78425.90005093576,786760.7841777843,-1.5944788233134339*1e7]
# zHDerivs_MMA = [1.,-0.6411496334238431,-11.243389142737914,-43.77408426449627,-106.78248444718885,11.278780437300327,2326.929754644706,19606.415106847817,113533.89599864538,464196.5937074214]
# xH_MMA = [xHDerivs_MMA[1], yHDerivs_MMA[1], zHDerivs_MMA[1]]
# dxH_MMA = [xHDerivs_MMA[2], yHDerivs_MMA[2], zHDerivs_MMA[2]]
# d2xH_MMA = [xHDerivs_MMA[3], yHDerivs_MMA[3], zHDerivs_MMA[3]]
# d3xH_MMA = [xHDerivs_MMA[4], yHDerivs_MMA[4], zHDerivs_MMA[4]]
# d4xH_MMA = [xHDerivs_MMA[5], yHDerivs_MMA[5], zHDerivs_MMA[5]]
# d5xH_MMA = [xHDerivs_MMA[6], yHDerivs_MMA[6], zHDerivs_MMA[6]]
# d6xH_MMA = [xHDerivs_MMA[7], yHDerivs_MMA[7], zHDerivs_MMA[7]]
# d7xH_MMA = [xHDerivs_MMA[8], yHDerivs_MMA[8], zHDerivs_MMA[8]]
# d8xH_MMA = [xHDerivs_MMA[9], yHDerivs_MMA[9], zHDerivs_MMA[9]]
# d9xH_MMA = [xHDerivs_MMA[10], yHDerivs_MMA[10], zHDerivs_MMA[10]]

# HarmonicCoordDerivs.compute_harmonic_derivs!(xBL, dxBL, d2xBL, d3xBL, d4xBL, d5xBL, d6xBL, d7xBL, d8xBL, d9xBL, xH, dxH, d2xH, d3xH, d4xH, d5xH, d6xH, d7xH, d8xH, d9xH, a)

# using Test
# @test xH ≈ xH_MMA
# @test dxH ≈ dxH_MMA
# @test d2xH ≈ d2xH_MMA
# @test d3xH ≈ d3xH_MMA
# @test d4xH ≈ d4xH_MMA
# @test d5xH ≈ d5xH_MMA
# @test d6xH ≈ d6xH_MMA
# @test d7xH ≈ d7xH_MMA
# @test d8xH ≈ d8xH_MMA
# @test d9xH ≈ d9xH_MMA
