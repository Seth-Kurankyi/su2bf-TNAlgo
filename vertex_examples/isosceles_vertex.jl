#include module(s)
include("../src/coherent_vertex.jl")
using .coherent_vertex


# use the following packages for the computation 
using JLD2, Memoize, BenchmarkTools

# face normal vectors for an equilateral tetrahedron
#[1,1,1,1]
n1= [[0.0, 0.0, 1.0], [0.0, 2sqrt(2)/3, -1/3], 
    [sqrt(2/3), -sqrt(2)/3, -1/3], 
    [-sqrt(2/3), -sqrt(2)/3, -1/3]];

#[1,1,1,2]
n2 = [[0.0, 0.0, 1.0],
 [0.0, 0.9860132971832694, 1/6],
 [0.9759000729485332, 0.14085904245475284, 1/6],
 [-0.4879500364742665, -0.563436169819011, -2/3]];


 bdyars = [1,1,1,2,1,1,2,1,2,2] 

function cohn_vertex(j) 
    return cohn_vertex2(j*bdyars,[n2,n2,n2,n2,n1])
end

println("---------- Initialize computation time ------------ ")
@time cohn_vertex(1.0);

println("----------- Start vertex computations ------ ")

spins = 0.0:20.0

isodataF = []
isodataC = []
@simd for i in spins 
    ss = @timed cohn_vertex(i)
    push!(isodataF,[i,ss.value,ss.time,ss.bytes] )
    print(" spin j = $i,"," ftime = ",ss.time)
    
    # ss = @timed cohn_vertex(i)
    # push!(isodataC,[i,ss.value,ss.time,ss.bytes] )
    # print(", ctime = ",ss.time)
    println(" ")

    # empty memoize functions to free up memory space 
    empty!(memoize_cache(wig3j))
    empty!(memoize_cache(wig6j))
    empty!(memoize_cache(wignermatrix))
end


# save data using JLD2

k1,k2=spins[1],spins[end]

@save "isovertex_$k1-to-$k2.jld2" isodataF isodataC
