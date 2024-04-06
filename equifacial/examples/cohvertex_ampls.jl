#include module(s)
include("../semi_coherent_vertices.jl")
using .semi_coherent_vertices


# use the following packages for the computation 
using JLD2, Memoize, BenchmarkTools

# face normal vectors for an equilateral tetrahedron
nn= [[0.0, 0.0, 1.0], [0.0, 2sqrt(2)/3, -1/3], 
    [sqrt(2/3), -sqrt(2)/3, -1/3], 
    [-sqrt(2/3), -sqrt(2)/3, -1/3]];


println("---------- Initialize computation time ------------ ")
@time cohn_vertex(1,nn);

println("----------- Start vertex computations ------ ")

spins = 0.0:0.5:100.0

eqdataF = []
eqdataC = []
@simd for i in spins 
    ss = @timed cohn_vertex(i,nn)
    push!(eqdataF,[i,ss.value,ss.time,ss.bytes] )
    print(" spin j = $i,"," ftime = ",ss.time," ,")
    
    ss = @timed cohn_vertex(i,nn)
    push!(eqdataC,[i,ss.value,ss.time,ss.bytes] )
    print(" ctime = ",ss.time)
    println(" ")

    # empty memoize functions to free up memory space 
    empty!(memoize_cache(wig3j))
    empty!(memoize_cache(wig6j))
    empty!(memoize_cache(wignermatrix))
end


# save data using JLD2

k1,k2=spins[1],spins[end]

@save "equivertex_$k1-to-$k2.jld2" eqdataF eqdataC
