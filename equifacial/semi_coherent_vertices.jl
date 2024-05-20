"""
Module `partial_coherent_vertex` compute the coherent vertex amplitude for equal spin assignments and normal vectors 
    The 'semi' refers to mixture of boundary and bulk edges or intertwiners 
"""

#change semi to partial
module semi_coherent_vertices

#include su2 functions from src 
include("../src/su2functions.jl")

using .su2functions

#export functions 
export  cohn_vertex,
        wig3j,
        wig6j,
        wignermatrix


"""
Compute the coherent vetex amplitude for equal spins j and normal vectors nv  
"""

function cohn_vertex(j,nv)
    # compute coherent 4j vector (with a phase)
    vc = vector_coherent4jPh(j*ones(4),nv)
    sol = 0.0im
    @simd for x in 3j:-1.0:0.0
        w6j = wig6j_matrix((0.0:2j,0.0:2j),j,x,j,j)
        w6j =  w6j .* vc 
        sol += (-1)^(2j)*(dim_j(x))* sum(transpose(w6j) .* (w6j^4))
    end
    return sol
end


#end module 
end
