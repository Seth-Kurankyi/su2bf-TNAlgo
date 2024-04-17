#include module(s)
include("../src/coherent_vertex.jl")
using .coherent_vertex


# use the following packages for the computation 
using JLD2, Memoize, BenchmarkTools

# face normal vectors for an equilateral tetrahedron
nn= [[0.0, 0.0, 1.0], [0.0, 2sqrt(2)/3, -1/3], 
    [sqrt(2/3), -sqrt(2)/3, -1/3], 
    [-sqrt(2/3), -sqrt(2)/3, -1/3]];


function cohn_vertex(j) 
    return cohn_vertex1(j*ones(10),[nn,nn,nn,nn,nn])
end

println("---------- Initialize computation time ------------ ")
@time cohn_vertex(1);

println("----------- Start vertex computations ------ ")

spins = 70.5:0.5:90.0

eqdataF = []
eqdataC = []
@simd for i in spins 
    ss = @timed cohn_vertex(i)
    push!(eqdataF,[i,ss.value,ss.time,ss.bytes] )
    print(" spin j = $i,"," ftime = ",ss.time)
    
    ss = @timed cohn_vertex(i)
    push!(eqdataC,[i,ss.value,ss.time,ss.bytes] )
    print(", ctime = ",ss.time)
    println(" ")

    # empty memoize functions to free up memory space 
    empty!(memoize_cache(wig3j))
    empty!(memoize_cache(wig6j))
    empty!(memoize_cache(wignermatrix))
end


# save data using JLD2

k1,k2=spins[1],spins[end]

@save "equivertex_$k1-to-$k2.jld2" eqdataF eqdataC
