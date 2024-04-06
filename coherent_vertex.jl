
"""
Compute the coherent vertex amplitude for generic spin assignments in terms ofmatrix contractions 
"""
# use functions from su2functions 
include("su2functions.jl")
using .su2functions


"""
Compute the coherent vetex amplitude for fixed virtual spin x: 
    the coherent 4j vectors are considered as inputs
"""
function cohnX_vertex(x,jays,vcs)
    #note the order of the spins 
    j12,j13,j14,j15,j23,j24,j25,j34,j35,j45 = jays

    #the groups of 4j intertwiner ranges must be consistent with the spin assignments 
    I1,I2,I3,I4,I5 = intw_range(j12,j13,j14,j15),intw_range(j12,j25,j24,j23),intw_range(j23,j13,j35,j34),
    intw_range(j34,j24,j14,j45),intw_range(j15,j25,j35,j45)
    
    # multiply components (of the first index) of the wigner 6j matrix by the coherent 4j vector (for the intertwiner index)
    w61 = wig6j_matrix([I1,I2],j25,x,j13,j12)
    w61 = w61 .* vcs[1]

    w62 = wig6j_matrix([I2,I3],j13,x,j24,j23)
    w62 = w62 .* vcs[2]

    w63 = wig6j_matrix([I3,I4],j24,x,j35,j34)
    w63 = w63 .* vcs[3]

    w64 = wig6j_matrix([I4,I5],j35,x,j14,j45)
    w64 = w64 .* vcs[4]

    w65 = wig6j_matrix([I5,I1],j14,x,j25,j15)
    w65 = w65 .* vcs[5]
    
    # compute trace of matrix multiplications for the intertwiners 
    return sum(transpose(w61) .* (w62*w63*w64*w65))
end




