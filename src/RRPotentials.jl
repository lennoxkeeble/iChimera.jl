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

function ε(i::Int, j::Int, k::Int)::Int
    return get(levi_civita_table, (i, j, k), 0)
end


function compute_RR_potentials!(Virr::MVector{3, Float64}, ∂Vrr_∂a::MVector{3, Float64}, ∂Virr_∂t::MVector{3, Float64}, ∂Virr_∂a::MMatrix{3, 3, Float64, 9}, x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, Mij5::MMatrix{3, 3, Float64, 9}, Mij6::MMatrix{3, 3, Float64, 9}, Mij7::MMatrix{3, 3, Float64, 9}, Mij8::MMatrix{3, 3, Float64, 9}, dxmMij5::MArray{Tuple{3, 3, 3}, Float64, 3, 27}, dxmMij6::MArray{Tuple{3, 3, 3}, Float64, 3, 27}, dxmMij7::MArray{Tuple{3, 3, 3}, Float64, 3, 27}, Mijk7::MArray{Tuple{3, 3, 3}, Float64, 3, 27}, Mijk8::MArray{Tuple{3, 3, 3}, Float64, 3, 27}, dxmMijk7::MArray{Tuple{3, 3, 3, 3}, Float64, 4, 81}, Skl5::MMatrix{3, 3, Float64, 9}, Skl6::MMatrix{3, 3, Float64, 9}, dxmSkl5::MArray{Tuple{3, 3, 3}, Float64, 3, 27}, n::Int64)

    Vrr = RRPotentials.Vrr(x, dx, Mij5, Mij7, Mijk7, n)
    ∂Vrr_∂t = RRPotentials.∂Vrr_∂t(x, dx, d2x, Mij5, Mij6, Mij7, Mij8, Mijk7, Mijk8, n)

    for i = 1:3
        ∂Vrr_∂a[i] = RRPotentials.∂Vrr_∂xm(i, x, dx, Mij5, dxmMij5, Mij7, dxmMij7, Mijk7, dxmMijk7, n)
        Virr[i] =  RRPotentials.Virr(i, x, Mij6, Skl5)
        ∂Virr_∂t[i] = RRPotentials.∂Virr_∂t(i, x, dx, Mij6, Mij7, Skl5, Skl6)
        for j = 1:3
            ∂Virr_∂a[i, j] = RRPotentials.∂Virr_∂xm(i, j, x, Mij6, dxmMij6, Skl5, dxmSkl5)
        end
    end 

    return Vrr, ∂Vrr_∂t
end

# function Vrr(x::MVector{3, Float64}, Mij5::MMatrix{3, 3, Float64, 9}, Mij7::MMatrix{3, 3, Float64, 9}, Mijk7::MArray{Tuple{3, 3, 3}, Float64, 3, 27})::Float64 
#     sum = 0.0
#     for i =1:3
#         for j =1:3
#             sum += (-1.0/5.0) * x[i] * x[j] * Mij5[i, j]
#             sum += (-1.0/70.0) * x[i] * x[j] * (x[1]^2+x[2]^2+x[3]^2) * Mij7[i, j]
#             for k =1:3
#                 sum += (x[i]*x[j]*x[k]*Mijk7[i,j,k])/189.
#             end
#         end
#     end
#     return sum
# end


# function ∂Vrr_∂t(x::MVector{3, Float64}, dx::MVector{3, Float64}, Mij5::MMatrix{3, 3, Float64, 9}, Mij6::MMatrix{3, 3, Float64, 9}, Mij7::MMatrix{3, 3, Float64, 9}, Mij8::MMatrix{3, 3, Float64, 9}, Mijk7::MArray{Tuple{3, 3, 3}, Float64, 3, 27}, Mijk8::MArray{Tuple{3, 3, 3}, Float64, 3, 27})::Float64 
#     sum = 0.0
#     for i =1:3
#         for j =1:3
#             sum += (-((x[j]*dx[i] + x[i]*dx[j])*Mij5[i,j]) - x[i]*x[j]*Mij6[i,j])/5.
#             sum += (-((x[j]*(2*x[i]*(x[1]*dx[1] + x[2]*dx[2] + x[3]*dx[3]) + (x[1]^2 + x[2]^2 + x[3]^2)*dx[i]) + (x[1]^2 + x[2]^2 + x[3]^2)*x[i]*dx[j])*Mij7[i,j]) - (x[1]^2 + x[2]^2 + x[3]^2)*x[i]*x[j]*Mij8[i,j])/70.
#             for k = 1:3
#                 sum += ((x[i]*x[k]*dx[j] + x[j]*(x[k]*dx[i] + x[i]*dx[k]))*Mijk7[i,j,k] + x[i]*x[j]*x[k]*Mijk8[i,j,k])/189.
#             end
#         end
#     end
#     return sum
# end

# function ∂Vrr_∂xm(m::Int64, x::MVector{3, Float64}, Mij5::MMatrix{3, 3, Float64, 9}, dxmMij5::MArray{Tuple{3, 3, 3}, Float64, 3, 27}, Mij7::MMatrix{3, 3, Float64, 9}, dxmMij7::MArray{Tuple{3, 3, 3}, Float64, 3, 27}, Mijk7::MArray{Tuple{3, 3, 3}, Float64, 3, 27}, dxmMijk7::MArray{Tuple{3, 3, 3, 3}, Float64, 4, 81})::Float64 
#     sum = 0.0
#     for i =1:3
#         for j =1:3
#             sum += (-(dxmMij5[i,j,m]*x[i]*x[j]) - Mij5[i,j]*(x[j]*δ(i,m) + x[i]*δ(j,m)))/5.
#             sum += (-(dxmMij7[i,j,m]*(x[1]^2 + x[2]^2 + x[3]^2)*x[i]*x[j]) - Mij7[i,j]*(2*x[i]*x[j]*(x[1]*δ(1,m) + x[2]*δ(2,m) + x[3]*δ(3,m)) + (x[1]^2 + x[2]^2 + x[3]^2)*(x[j]*δ(i,m) + x[i]*δ(j,m))))/70.
#             for k = 1:3
#                 sum += (dxmMijk7[i,j,k,m]*x[i]*x[j]*x[k] + Mijk7[i,j,k]*(x[i]*x[k]*δ(j,m) + x[j]*(x[k]*δ(i,m) + x[i]*δ(k,m))))/189. 
#             end
#         end
#     end
#     return sum
# end


function Vrr(x::MVector{3, Float64}, dx::MVector{3, Float64}, Mij5::MMatrix{3, 3, Float64, 9}, Mij7::MMatrix{3, 3, Float64, 9}, Mijk7::MArray{Tuple{3, 3, 3}, Float64, 3, 27}, n::Int64)::Float64 
    sum = 0.0
    for i =1:3
        for j =1:3
            sum += (-1.0/5.0) * x[i] * x[j] * Mij5[i, j]
            sum += (dx[1]^2 + dx[2]^2 + dx[3]^2)^n * (-1.0/70.0) * x[i] * x[j] * (x[1]^2+x[2]^2+x[3]^2) * Mij7[i, j]
            for k =1:3
                sum += (dx[1]^2 + dx[2]^2 + dx[3]^2)^n * (x[i]*x[j]*x[k]*Mijk7[i,j,k])/189.
            end
        end
    end
    return sum
end


function ∂Vrr_∂t(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, Mij5::MMatrix{3, 3, Float64, 9}, Mij6::MMatrix{3, 3, Float64, 9}, Mij7::MMatrix{3, 3, Float64, 9}, Mij8::MMatrix{3, 3, Float64, 9}, Mijk7::MArray{Tuple{3, 3, 3}, Float64, 3, 27}, Mijk8::MArray{Tuple{3, 3, 3}, Float64, 3, 27}, n::Int64)::Float64 
    sum = 0.0
    for i =1:3
        for j =1:3
            sum += (-((x[j]*dx[i] + x[i]*dx[j])*Mij5[i,j]) - x[i]*x[j]*Mij6[i,j])/5.
            sum += ((dx[1]^2 + dx[2]^2 + dx[3]^2)^n*(-2*x[i]*x[j]*(x[1]*dx[1] + x[2]*dx[2] + x[3]*dx[3])*Mij7[i,j] - (x[1]^2 + x[2]^2 + x[3]^2)*x[j]*dx[i]*Mij7[i,j] - (x[1]^2 + x[2]^2 + x[3]^2)*x[i]*dx[j]*Mij7[i,j] - (2*n*(x[1]^2 + x[2]^2 + x[3]^2)*x[i]*x[j]*(dx[1]*d2x[1] + dx[2]*d2x[2] + dx[3]*d2x[3])*Mij7[i,j])/(dx[1]^2 + dx[2]^2 + dx[3]^2) - (x[1]^2 + x[2]^2 + x[3]^2)*x[i]*x[j]*Mij8[i,j]))/70.
            for k = 1:3
                sum += ((dx[1]^2 + dx[2]^2 + dx[3]^2)^n*(x[j]*x[k]*dx[i]*Mijk7[i,j,k] + x[i]*x[k]*dx[j]*Mijk7[i,j,k] + x[i]*x[j]*dx[k]*Mijk7[i,j,k] + (2*n*x[i]*x[j]*x[k]*(dx[1]*d2x[1] + dx[2]*d2x[2] + dx[3]*d2x[3])*Mijk7[i,j,k])/(dx[1]^2 + dx[2]^2 + dx[3]^2) + x[i]*x[j]*x[k]*Mijk8[i,j,k]))/189.
            end
        end
    end
    return sum
end

function ∂Vrr_∂xm(m::Int64, x::MVector{3, Float64}, dx::MVector{3, Float64}, Mij5::MMatrix{3, 3, Float64, 9}, dxmMij5::MArray{Tuple{3, 3, 3}, Float64, 3, 27}, Mij7::MMatrix{3, 3, Float64, 9}, dxmMij7::MArray{Tuple{3, 3, 3}, Float64, 3, 27}, Mijk7::MArray{Tuple{3, 3, 3}, Float64, 3, 27}, dxmMijk7::MArray{Tuple{3, 3, 3, 3}, Float64, 4, 81}, n::Int64)::Float64 
    sum = 0.0
    for i =1:3
        for j =1:3
            sum += (-(dxmMij5[i,j,m]*x[i]*x[j]) - Mij5[i,j]*(x[j]*δ(i,m) + x[i]*δ(j,m)))/5.
            sum += -0.014285714285714285*((dxmMij7[i,j,m]*(x[1]^2 + x[2]^2 + x[3]^2)*x[i]*x[j] + Mij7[i,j]*(x[j]*(2*x[i]*(x[1]*δ(1,m) + x[2]*δ(2,m) + x[3]*δ(3,m)) + (x[1]^2 + x[2]^2 + x[3]^2)*δ(i,m)) + (x[1]^2 + x[2]^2 + x[3]^2)*x[i]*δ(j,m)))*(dx[1]^2 + dx[2]^2 + dx[3]^2)^n)
            for k = 1:3
                sum += ((dxmMijk7[i,j,k,m]*x[i]*x[j]*x[k] + Mijk7[i,j,k]*(x[i]*x[k]*δ(j,m) + x[j]*(x[k]*δ(i,m) + x[i]*δ(k,m))))*(dx[1]^2 + dx[2]^2 + dx[3]^2)^n)/189.
            end
        end
    end
    return sum
end

function Virr(i::Int64, x::MVector{3, Float64}, Mij6::MMatrix{3, 3, Float64, 9}, Skl5::MMatrix{3, 3, Float64, 9})::Float64 
    sum = 0.0
    for j =1:3
        for k =1:3
            sum += ((x[i]*x[j]*x[k] - ((x[1]^2 + x[2]^2 + x[3]^2)*(x[k]*δ(i,j) + x[j]*δ(i,k) + x[i]*δ(j,k)))/5.)*Mij6[j,k])/21.
            for l = 1:3
                sum += (-4*x[j]*x[l]*ε(i,j,k)*Skl5[k,l])/45.
            end
        end
    end
    return sum
end

function ∂Virr_∂t(i::Int64, x::MVector{3, Float64}, dx::MVector{3, Float64}, Mij6::MMatrix{3, 3, Float64, 9}, Mij7::MMatrix{3, 3, Float64, 9}, Skl5::MMatrix{3, 3, Float64, 9}, Skl6::MMatrix{3, 3, Float64, 9})::Float64 
    sum = 0.0
    for j =1:3
        for k =1:3
            sum += (((-2*(x[k]*δ(i,j) + x[j]*δ(i,k) + x[i]*δ(j,k))*(x[1]*dx[1] + x[2]*dx[2] + x[3]*dx[3]))/5. + x[j]*x[k]*dx[i] + x[i]*x[k]*dx[j] + x[i]*x[j]*dx[k] - ((x[1]^2 + x[2]^2 + x[3]^2)*(δ(j,k)*dx[i] + δ(i,k)*dx[j] + δ(i,j)*dx[k]))/5.)*Mij6[j,k] + (x[i]*x[j]*x[k] - ((x[1]^2 + x[2]^2 + x[3]^2)*(x[k]*δ(i,j) + x[j]*δ(i,k) + x[i]*δ(j,k)))/5.)*Mij7[j,k])/21.
            for l = 1:3
                sum += (-4*ε(i,j,k)*((x[l]*dx[j] + x[j]*dx[l])*Skl5[k,l] + x[j]*x[l]*Skl6[k,l]))/45.
            end
        end
    end
    return sum
end

function ∂Virr_∂xm(i::Int64, m::Int64, x::MVector{3, Float64}, Mij6::MMatrix{3, 3, Float64, 9}, dxmMij6::MArray{Tuple{3, 3, 3}, Float64, 3, 27}, Skl5::MMatrix{3, 3, Float64, 9}, dxmSkl5::MArray{Tuple{3, 3, 3}, Float64, 3, 27})::Float64 
    sum = 0.0
    for j =1:3
        for k =1:3
            sum += (dxmMij6[j,k,m]*(x[i]*x[j]*x[k] - ((x[1]^2 + x[2]^2 + x[3]^2)*(x[k]*δ(i,j) + x[j]*δ(i,k) + x[i]*δ(j,k)))/5.) + Mij6[j,k]*(x[j]*x[k]*δ(i,m) - (2*(x[1]*δ(1,m) + x[2]*δ(2,m) + x[3]*δ(3,m))*(x[k]*δ(i,j) + x[j]*δ(i,k) + x[i]*δ(j,k)))/5. + x[i]*(x[k]*δ(j,m) + x[j]*δ(k,m)) - ((x[1]^2 + x[2]^2 + x[3]^2)*(δ(i,m)*δ(j,k) + δ(i,k)*δ(j,m) + δ(i,j)*δ(k,m)))/5.))/21.
            for l = 1:3
                sum += (-4*(dxmSkl5[k,l,m]*x[j]*x[l] + Skl5[k,l]*(x[l]*δ(j,m) + x[j]*δ(l,m)))*ε(i,j,k))/45.
            end
        end
    end
    return sum
end

end

# using StaticArrays
# using Random
# x = @MVector rand(3)
# dx = @MVector rand(3)
# Mij5 = @MArray rand(3, 3)
# Mij6 = @MArray rand(3, 3)
# Mij7 = @MArray rand(3, 3)
# Mij8 = @MArray rand(3, 3)
# dxmMij5 = @MArray rand(3, 3, 3)
# dxmMij6 = @MArray rand(3, 3, 3)
# dxmMij7 = @MArray rand(3, 3, 3)

# Mijk7 = @MArray rand(3, 3, 3)
# Mijk8 = @MArray rand(3, 3, 3)
# dxmMijk7 = @MArray rand(3, 3, 3, 3)

# Skl5 = @MArray rand(3, 3)
# Skl6 = @MArray rand(3, 3)
# dxmSkl5 = @MArray rand(3, 3, 3)

# Virr = @MVector zeros(3)
# ∂Vrr_∂a = @MVector zeros(3)
# ∂Virr_∂t = @MVector zeros(3)
# ∂Virr_∂a = @MArray zeros(3, 3)

# RRPotentials.compute_RR_potentials!(Virr, ∂Vrr_∂a, ∂Virr_∂t, ∂Virr_∂a, x, dx, Mij5, Mij6, Mij7, Mij8, dxmMij5, dxmMij6, dxmMij7, Mijk7, Mijk8, dxmMijk7, Skl5, Skl6, dxmSkl5)