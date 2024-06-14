"""
Module `coherent_vertex` compute the coherent vertex amplitude for generic spin assignments and normal vectors in terms of matrix contractions
"""

module coherent_vertex

# use functions from su2functions 
include("su2functions.jl")
using .su2functions


#export the following 

export  wig3j, 
        wig6j, 
        wignermatrix, 
        cohnX_vertex,
        cohn_vertex1,
        cohn_vertex2


"""
Compute the coherent vetex amplitude for fixed virtual spin x: 
    the coherent {4j} vectors are considered as inputs
"""
function cohnX_vertex(x,jays,vcs)
    #note the order of the spins, 1,2,3,4,5 labels the edges (tetrahedron)
    # the labels are equivalent to the bar notation where 1,2,3,4,5 label vertices.. 
    j12, j13, j14, j15, j23, j24, j25, j34, j35, j45 = jays

    #the groups of 4j intertwiner ranges must be consistent with the spin assignments .. 
    #clockwise order for all intertwiners according to graph of vertex
    I1,I2,I3,I4,I5 = intw_range(j12 ,j13 ,j14 ,j15), intw_range(j23, j24, j25, j12), intw_range(j34,j35, j13, j23),
    intw_range(j45, j14, j24, j34), intw_range(j15, j25, j35, j45)

    #coherent 6j matrices
    f6j12 = coherent_wig6j_matrix(vc[1],(I1,I2),j25,x,j13,j12)

    f6j23 = coherent_wig6j_matrix(vc[2],(I2,I3),j13,x,j24,j23)

    f6j34 = coherent_wig6j_matrix(vc[3],(I3,I4),j24,x,j35,j34)

    f6j45 = coherent_wig6j_matrix(vc[4],(I4,I5),j35,x,j14,j45)

    f6j51 = coherent_wig6j_matrix(vc[5],(I5,I1),j14,x,j25,j15)

    # compute trace of matrix multiplications for the intertwiners 
    return sum( transpose(f6j612) .* (f6j23 * f6j34 * f6j45 * f6j51) )
end


# function cohnX_vertex(x,jays,vcs)
#     #note the order of the spins, 1,2,3,4,5 labels the edges (tetrahedron)
#     # the labels are equivalent to the bar notation where 1,2,3,4,5 label vertices.. 
#     j12, j13, j14, j15, j23, j24, j25, j34, j35, j45 = jays

#     #the groups of 4j intertwiner ranges must be consistent with the spin assignments .. 
#     #clockwise order for all intertwiners according to graph of vertex
#     I1,I2,I3,I4,I5 = intw_range(j12 ,j13 ,j14 ,j15), intw_range(j23, j24, j25, j12), intw_range(j34,j35, j13, j23),
#     intw_range(j45, j14, j24, j34), intw_range(j15, j25, j35, j45)
    
#     # multiply components (of the first index) of the wigner 6j matrix by the coherent 4j vector (for the intertwiner index)
#     w612 = wig6j_matrix((I1,I2),j25,x,j13,j12)
#     w612 = w612 .* vcs[1]

#     w623 = wig6j_matrix((I2,I3),j13,x,j24,j23)
#     w623 = w623 .* vcs[2]

#     w634 = wig6j_matrix((I3,I4),j24,x,j35,j34)
#     w634 = w634 .* vcs[3]

#     w645 = wig6j_matrix((I4,I5),j35,x,j14,j45)
#     w645 = w645 .* vcs[4]

#     w651 = wig6j_matrix((I5,I1),j14,x,j25,j15)
#     w651 = w651 .* vcs[5]
    
#     # compute trace of matrix multiplications for the intertwiners 
#     return sum( transpose(w612) .* (w623 * w634 * w645 * w651) )
# end


"""
Compute the coherent vertex amplitude for given the spins and normal vectors as inputs
"""
function cohn_vertex1(jays,nvs)
    j12,j13,j14,j15,j23,j24,j25,j34,j35,j45 = jays
    nv1,nv2,nv3,nv4,nv5 = nvs
    
    # !!!!!!!   IMPORTANT   !!!!!!!
    """ The order of the spins are crucial for the matching the coherent vectors and the input normal vectors: 
        Note that order of the spins should match the intertwiners and also match the order of the normal vectors (inputs)
    """
    # The order of the spins are the same as the order of the unit normal vectors as inputs
    jjs = (j12,j13,j14,j15),(j23,j24,j25,j12),(j34,j35,j13,j23),(j45,j14,j24,j34),(j15,j25,j35,j45)

        
    # intertwiner ranges 
    #I1,I2,I3,I4,I5 = intw_range.(jjs)

    # compute the coherent {4j} vector for all 5 boundary edges 
    vcs = Dict()
    for i in 1:5 
        vcs[i] = vector_coherent4jPh(jjs[i],nvs[i])
    end
    
    #vcs = vc1,vc2,vc3,vc4,vc5
    
    sol = 0.0im 
    # range of values for the virtual spin 
    x_range = min( min(j12+j13,j14+j15)+j23, min(j24+j23,j12+j25)+j13, min(j13+j23,j34+j35)+j24,
                min(j34+j24,j14+j45)+j35 , min(j15+j25,j35+j45)+j14 )
    
    @simd for x in x_range:-1:0
        sol += (-1.0+0.0im)^(sum(jays))*(dim_j(x))*cohnX_vertex(x,jays,vcs)
    end
    return sol
end

"""
Coherent vertex amplitude for given the spins and normal vectors as inputs with different ordering 
"""
function cohn_vertex2(jays,nvs)
    j12,j13,j14,j15,j23,j24,j25,j34,j35,j45 = jays
    nv1,nv2,nv3,nv4,nv5 = nvs
    
    # !!!!!!!   IMPORTANT   !!!!!!!
    """ The order of the spins are crucial for the matching the coherent vectors and the input normal vectors: 
        Note that order of the spins should match the intertwiners and also match the order of the normal vectors (inputs)
    """
    # The order of the spins are the same as the order of the unit normal vectors as inputs
    jjs = (j12,j13,j14,j15),(j24,j23,j12,j25),(j13,j23,j34,j35),(j34,j24,j14,j45),(j15,j25,j35,j45)

        
    # intertwiner ranges 
    #I1,I2,I3,I4,I5 = intw_range.(jjs)

    # compute the coherent {4j} vector for all 5 boundary edges 
    vcs = Dict()
    for i in 1:5 
        vcs[i] = vector_coherent4jPh(jjs[i],nvs[i])
    end
    
    #vcs = vc1,vc2,vc3,vc4,vc5
    
    sol = 0.0im 
    # range of values for the virtual spin 
    x_range = min( min(j12+j13,j14+j15)+j23, min(j24+j23,j12+j25)+j13, min(j13+j23,j34+j35)+j24,
                min(j34+j24,j14+j45)+j35 , min(j15+j25,j35+j45)+j14 )
    
    @simd for x in x_range:-1:0
        sol += (-1.0+0.0im)^(sum(jays))*(dim_j(x))*cohnX_vertex(x,jays,vcs)
    end
    return sol
end


#end module 
end
