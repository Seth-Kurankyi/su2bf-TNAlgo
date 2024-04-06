"""
Module `su2functions` contains utility functions for working with SU(2) representation theory functions including Wigner symbols.
"""

module su2functions

using LinearAlgebra, Combinatorics, HalfIntegers, WignerSymbols, Memoize

# Exported functions
export  is_ang, 
        is_pair, 
        intw_range, 
        dim_j, 
        wig3j, 
        wig6j, 
        wignermatrix, 
        coherent4j, 
        vector_coherent4j, 
        vector_coherent4jPh
        wig6j_matrix, 

#======== Utility functions for SU(2) angular momentum and Wigner symbols ============#

"""
Check if a spin and magnetic indices (j,m) for an angular momentum pair
"""
function is_ang(j, m)::Bool
    return abs(m) ≤ j && ishalfinteger(j) && isinteger(j - m) && isinteger(j + m)
end

""" 
Check intergness of  pair of spins (j1,j2) 
"""
function is_pair(j1, j2)::Bool
    return ishalfinteger(j1) && isinteger(j1 - j2) && isinteger(j1 + j2)
end

"""
Compute the intertwiner range for four spins representations 
"""
function intw_range(j1, j2, j3, j4)
    return max(abs(j1 - j2), abs(j3 - j4)) : min(j1 + j2, j3 + j4)
end

# Convert a vector of spins to an intertwiner range
function intw_range(js)
    return intw_range(js[1], js[2], js[3], js[4])
end

"""
Compute the dimension of spin j
"""
dim_j(j)::Float64 = twice(j) + 1.0

#========= Functions for Wigner matrices of a unit normal vector n ==============#

"""
Compute the S^2 angles (ϕ,θ) from a 3D (unit) normal vector
"""
function s2angles(n)::Vector{Float64}
    if !isapprox(norm(n),1.0 ) # Check if normal vector is normalized 
        n = normalize(n) # normalize
    end
    return [atan(n[2], n[1]), acos(n[3])] # Conditions for n[1] == 0 could be added
end

"""
Compute components of Wigner matrix as a function of the S^2 angles (ϕ,θ) in the |+> or highest spin weight sector
"""
function wignermatrix_ang(j, m, phi, theta)::ComplexF64
    sin_th2, cos_th2 = sin(theta/2), cos(theta/2)
    # sin(θ/2) =0  or cos(θ/2) = 0 
    if isapprox(sin_th2,0.0) # Check if sin(θ/2) ≈ 0.0
        return j == m ? 1.0 : 0.0
    elseif isapprox(cos_th2,0.0) # Check if cos(θ/2) ≈ 0.0
        return j == -m ? 1.0 : 0.0
    elseif m == 0.0 # for m=0 binomial function simplies in terms of catalan numbers 
        return sqrt((j + 1) * catalannum(Int(j))) * (cos_th2 * sin_th2 * exp(-im * phi)) ^ j
    elseif m == 1.0 || m == -1.0 # for m=1 binomial function simplies in terms of catalan numbers 
        return sqrt(j * catalannum(Int(j))) * (cos_th2) ^ (j + m) * (sin_th2 * exp(-im * phi)) ^ (j - m)
    else
        # For j > 33 use BigInt to avoid overflows in the binomial function
        j <= 33 ? two_j = Int(2 * j) : two_j = big(Int(2 * j))
        return sqrt(binomial(two_j, Int(j + m))) * (cos_th2) ^ (j + m) * (sin_th2 * exp(-im * phi)) ^ (j - m)
    end
end

"""
Compute components of Wigner matrix as a function of 3D unit normal vecto in the |+> 
"""
function wignermatrix_nvec(j, m, nv)::ComplexF64
    return wignermatrix_ang(j,m,atan(nv[2],nv[1]),acos(nv[3]))
end

"""
Memoized version of wignermatrix_nvec function
"""
@memoize function wignermatrix(j, m, nv)::ComplexF64
   return wignermatrix_nvec(j, m, nv)
end

# Functions for Wigner 3j symbols

"""
Memoized version of Wigner 3j symbol
"""
@memoize function wig3j(j1, j2, j3, m1, m2, m3 = -m1 - m2)::Float64
    return wigner3j(j1, j2, j3, m1, m2, m3)
end

# Functions for Wigner 4j symbols: 4j intertwiner at channel j12 - for all incoming edges

"""
Compute the Wigner 4j symbol
"""
function wigner4j(j12, j1, j2, j3, j4, m1, m2, m3, m4 = -m1 - m2 - m3)::Float64
    return is_ang(j12,m1+m2) ? (-1)^(j12+m1+m2)*wig3j(j1,j2,j12,m1,m2,-m1-m2)*wig3j(j12,j3,j4,m1+m2,m3,m4) : 0.0 
end

# Wigner 4j symbols in coherent basis

"""
Compute the coherent 4j intertwiner (channel iota): as a function of spins and normal vectors
"""
function coherent4j(iota, jay, nvecs)::ComplexF64
    r::ComplexF64 = 0.0
    j1, j2, j3, j4 = jay
    n1, n2, n3, n4 = nvecs
    for m1 in -j1:j1
        wm1 = wignermatrix(j1, m1, n1) 
        if wm1 != 0.0
            for m2 in -j2:j2
                if is_ang(iota, m1 + m2) 
                    wm2 = wignermatrix(j2, m2, n2)
                    if wm2 != 0.0
                        for m3 in max(-j3, -m1 - m2 - j4):min(j3, j4 - m1 - m2)
                            wm34 = wignermatrix(j3, m3, n3) * wignermatrix(j4, -(m1 + m2 + m3), n4)
                            w4j = wigner4j(iota, j1, j2, j3, j4, m1, m2, m3)
                            if wm34 != 0.0 && w4j != 0.0
                                r += wm1 * wm2 * wm34 * w4j
                            end
                        end
                    end
                end
            end
        end
    end
    return r
end

"""
Compute a vector of coherent 4j symbol for all possible values of iota
"""
function vector_coherent4j(js, nvecs)::Vector{ComplexF64}
    #js is a vector of 4 spins 
    itr = intw_range(js)
    ls = length(itr)
    #create a vector for the length of range of the intertwiner
    sol = Vector{ComplexF64}(undef, ls)
    @simd for i in 1:ls
        sol[i] = dim_j(itr[i])*coherent4j(itr[i],js,nvecs)
    end
    return sol
end

# Compute vector of coherent 4j amplitudes with phase
function vector_coherent4jPh(js, nvecs)::Vector{ComplexF64}
    #js is a vector of 4 spins 
    itr = intw_range(js)
    ls = length(itr)
    #create a vector for the length of range of the intertwiner
    sol = Vector{ComplexF64}(undef, ls)
    @simd for i in 1:ls
        sol[i] = (-1.0+0.0im)^(itr[i])*dim_j(itr[i])*coherent4j(itr[i],js,nvecs)
    end
    return sol
end

# Functions for Wigner 6j symbols

"""
Memoized version of Wigner 6j symbol
"""
@memoize function wig6j(j1, j2, j3, j4, j5, j6)::Float64
    return wigner6j(j1, j2, j3, j4, j5, j6)
end

"""
Compute Wigner 6j symbol (i1, ja, jb ; i2, jb, jc) as a matrix in the indices i1, i2   
"""
function wig6j_matrix(IRs, ja, x, jb, jc)
    I1,I2 = IRs # range of values for indices i1, i2 (intertwiners)
    l1,l2 = length(I1), length(I2)
    sol = Matrix{Float64}(undef, l1, l2) 
    if is_pair(x,ja) && is_pair(x,jb) && is_pair(ja,jc) && is_pair(jb,jc)
        if ja == jb && I1 == I2 # the matrix is symmetric
            ir1 = intw_range(x,ja,jb,jc)
            if !isempty(ir1)
                for i1 in 1:l1
                    if I1[i1] in ir1
                        @simd for i2 in i1:l1
                            sol[i2,i1]  = sol[i1,i2] = wig6j(I1[i1],ja,x,I1[i2],jb,jc)
                        end
                    else
                        fill!(view(sol, i1, :), 0.0)
                        fill!(view(sol, :, i1), 0.0)
                    end
                end
            else
                fill!(view(sol, :, :), 0.0)
            end
        else
            ir1 = intw_range(x,ja,jb,jc)
            ir2 = intw_range(x,jb,ja,jc)
             if !isempty(ir1) || !isempty(ir2)
                for i1 in 1:l1
                    if I1[i1] in ir1
                        @simd for i2 in 1:l2
                            if I2[i2] in ir2
                                sol[i1,i2] = wig6j(I1[i1],ja,x,I2[i2],jb,jc)
                            else
                                sol[i1,i2] = 0.0
                            end
                        end
                    else 
                        fill!(view(sol, i1, :), 0.0)
                    end
                end
            else
                fill!(view(sol, :, :), 0.0)
            end
        end
    else
        fill!(view(sol, :, :), 0.0)
    end
    return sol
end

end # module