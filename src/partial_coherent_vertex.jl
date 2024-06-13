"""
Module `partial_coherent_vertex` compute the coherent vertex amplitude for 'boundary vertices' using matrix contractions.
"""

module partial_coherent_vertex

# use functions from su2functions 
include("su2functions.jl")
using .su2functions

#export the following functions 
#export  

"""
Compute the partial coherent vertex amplitude for fixed virtual spin x: 
    the coherent {4j} vectors are considered as inputs
"""
function partial_cohnX_vertex1(x,jays,vcs)
    #note the order of the spins, 1,2,3,4,5 labels the edges (tetrahedron)
    # the labels are equivalent to the bar notation where 1,2,3,4,5 label vertices.. 
    j12, j13, j14, j15, j23, j24, j25, j34, j35, j45 = jays

    #the groups of 4j intertwiner ranges must be consistent with the spin assignments .. 
    #clockwise order for all intertwiners according to graph of vertex
    I1,I2,I3,I4,I5 = intw_range(j12 ,j13 ,j14 ,j15), intw_range(j23, j24, j25, j12), intw_range(j34,j35, j13, j23),
    intw_range(j45, j14, j24, j34), intw_range(j15, j25, j35, j45)
    
    # multiply components (of the first index) of the wigner 6j matrix by the coherent 4j vector (for the intertwiner index)
    w612 = wig6j_matrix((I1,I2),j25,x,j13,j12)

    w623 = wig6j_matrix((I2,I3),j13,x,j24,j23)
    w623 = w623 .* vcs[1]

    w634 = wig6j_matrix((I3,I4),j24,x,j35,j34)
    w634 = w634 .* vcs[2]

    w645 = wig6j_matrix((I4,I5),j35,x,j14,j45)
    w645 = w645 .* vcs[3]

    w651 = wig6j_matrix((I5,I1),j14,x,j25,j15)
    w651 = w651 .* vcs[4]
    
    # compute sum over intertwiners as matrix multiplications 
    return transpose(sum( (w612) .* transpose(w623 * w634 * w645 * w651) , dims = 1))
end


function partial_cohnX_vertex2(x,jays,vcs)
    #note the order of the spins, 1,2,3,4,5 labels the edges (tetrahedron)
    # the labels are equivalent to the bar notation where 1,2,3,4,5 label vertices.. 
    j12, j13, j14, j15, j23, j24, j25, j34, j35, j45 = jays

    #the groups of 4j intertwiner ranges must be consistent with the spin assignments .. 
    #clockwise order for all intertwiners according to graph of vertex
    I1,I2,I3,I4,I5 = intw_range(j12 ,j13 ,j14 ,j15), intw_range(j23, j24, j25, j12), intw_range(j34,j35, j13, j23),
    intw_range(j45, j14, j24, j34), intw_range(j15, j25, j35, j45)
    
    # multiply components (of the first index) of the wigner 6j matrix by the coherent 4j vector (for the intertwiner index)
    w612 = wig6j_matrix((I1,I2),j25,x,j13,j12)

    w623 = wig6j_matrix((I2,I3),j13,x,j24,j23)

    w634 = wig6j_matrix((I3,I4),j24,x,j35,j34)
    w634 = w634 .* vcs[1]

    w645 = wig6j_matrix((I4,I5),j35,x,j14,j45)
    w645 = w645 .* vcs[2]

    w651 = wig6j_matrix((I5,I1),j14,x,j25,j15)
    w651 = w651 .* vcs[3]
    
    # compute sum over intertwiners as matrix multiplications 
    # returns a matrix
    return (w612) .* transpose(w623 * w634 * w645 * w651) 
end

end