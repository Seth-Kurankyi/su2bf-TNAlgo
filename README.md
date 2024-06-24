This repository contains functions for computing spin foam amplitudes for SU(2) BF theory. The algorithms here improve the computation of spin foam amplitudes via tensor contractions of low-valence tensors. They are solely written in **Julia** language. 

Spin foam amplitudes are typically associated with a 2-complex which is dual to a triangulation. 


### Dependencies
* The SU(2) invariant functions rely on [HalfIntegers.jl](https://github.com/sostock/HalfIntegers.jl) and [WignerSymbols.jl](https://github.com/Jutho/WignerSymbols.jl).
* We have also used [Memoize.jl](https://github.com/JuliaCollections/Memoize.jl) to optimize computations of certain functions through caching.


### Vertex amplitudes:
Coherent vertex amplitudes are computed using matrix contractions. The input are 10 spins $j_{ab}$ and a set of unit normal vectors $\bf n$ for each edge.
Here's a sample code to compute the coherent vertex amplitude with equal spins `ones(10)` and normal vectors `nn1,nn2,nn3,nn4,nn5`

`>> coherent_vertex1(ones(10),[nn1,nn2,nn3,nn4,nn5])`




