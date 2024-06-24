This repository contains functions for computing spin foam amplitudes for SU(2) BF theory. The algorithms here improve the computation of spin foam amplitudes via tensor contractions of low-valence tensors. They are solely written in **Julia** language. 

<!-- Spin foam amplitude is typically associated with a 2-complex which is dual to a triangulations. Tbe boundary data for the coherent amplitudes are spin labels and normal vectors. -->


### Dependencies
* The SU(2) invariant functions rely on [HalfIntegers.jl](https://github.com/sostock/HalfIntegers.jl) and [WignerSymbols.jl](https://github.com/Jutho/WignerSymbols.jl).
* We have also used [Memoize.jl](https://github.com/JuliaCollections/Memoize.jl) to optimize computations of certain functions through caching.


### Vertex amplitudes:
Coherent vertex amplitudes are computed using matrix contractions. The input are 10 spins $j_{ab}$ and a set of unit normal vectors $\bf n$ for each edge. It outputs a complex number. The arrangements of the normal vectors is very important. See the functions for more details. 

To compute the coherent vertex amplitude with equal spins `ones(10)` and normal vectors `nn1,nn2,nn3,nn4,nn5` do the following:

`>> coherent_vertex1(ones(10),[nn1,nn2,nn3,nn4,nn5])`


### Partial-Coherent amplitudes:
Partial coherent vertex amplitudes are associated with vertices with both bulk and boundary edges. They are useful for computing amplitudes for complexes with multiple vertices. 
The inputs are the 10 spins and the normal vectors for boundary edges. It ouputs a $n$-valent tensor if there are $n$ bulk edges. 

For example, the partial-coherent vertex amplitude with normals vectors `nn1,nn2,nn3` for boundary edges outputs a matrix with the bulk intertwiners as indices. 

`>> partial_cohn_vertex2(ones(10),nn1,nn2,nn3)`





