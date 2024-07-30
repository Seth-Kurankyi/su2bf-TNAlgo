#include module(s)
include("../src/coherent_vertex.jl")
include("../src/su2functions.jl")
using .coherent_vertex
using .su2functions


# use the following packages for the computation 

# face normal vectors for a degenerate tetrahedron
nv= [[0.0, 0.0, 1.0], [0.0, 0.0, -1.0], [0.0, 0.0, 1.0], [0.0, 0.0, -1.0]];

function vertex_degenerate(j) 
    return cohn_vertex1(j*ones(10),[nv,nv,nv,nv,nv])
end

println(vertex_degenerate(1))