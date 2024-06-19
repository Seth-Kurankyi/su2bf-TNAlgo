#add test soon 

include("../src/coherent_vertex.jl")

using .coherent_vertex

# face normal vectors for an equilateral tetrahedron
nn= [[0.0, 0.0, 1.0], [0.0, 2sqrt(2)/3, -1/3], 
    [sqrt(2/3), -sqrt(2)/3, -1/3], 
    [-sqrt(2/3), -sqrt(2)/3, -1/3]];



@time cohn_vertex1(ones(10),[nn,nn,nn,nn,nn])

