# SU(2) BF Spin Foam Amplitudes in Julia

This repository contains functions for computing spin foam amplitudes for SU(2) BF theory. The algorithms here improve the computation of spin foam amplitudes via tensor contractions of low-valence tensors. They are solely written in **Julia** language.

### Dependencies
* The SU(2) invariant functions rely on [HalfIntegers.jl](https://github.com/sostock/HalfIntegers.jl) and [WignerSymbols.jl](https://github.com/Jutho/WignerSymbols.jl).
* We have also used [Memoize.jl](https://github.com/JuliaCollections/Memoize.jl) to optimize computations of certain functions through caching.

### Vertex Amplitudes
Coherent vertex amplitudes are computed using matrix contractions. The input consists of ten spins \( j_{ab} \) and a set of unit normal vectors \( \{ \mathbf{n} \} \) for each edge. The function outputs a complex number. The arrangement of the normal vectors is crucial for the computation. See the functions for more details.

To compute the coherent vertex amplitude with spins `(j12, j13, j14, j15, j23, j24, j25, j34, j35, j45)` and set of normal vectors `(nn1, nn2, nn3, nn4, nn5)`, use the following command:

```julia
coherent_vertex1([j12, j13, j14, j15, j23, j24, j25, j34, j35, j45], [nn1, nn2, nn3, nn4, nn5])
```

### Partial-Coherent Amplitudes
Partial coherent vertex amplitudes are associated with vertices that have both bulk and boundary edges. They are useful for computing amplitudes for complexes with multiple vertices. The inputs are the 10 spins and the normal vectors for boundary edges. It outputs a \( k \)-valent tensor if there are \( k \) bulk edges.

For example, to compute the partial-coherent vertex amplitude with normal vectors `(nn1, nn2, nn3)` for boundary edges, which outputs a matrix with the bulk intertwiners as indices, use:

```julia
partial_cohn_vertex2([j12, j13, j14, j15, j23, j24, j25, j34, j35, j45], [nn1, nn2, nn3])
```

### Installation
To install the necessary dependencies and clone this repository, you can use the following commands in your Julia REPL:

```julia
import Pkg
Pkg.add("HalfIntegers")
Pkg.add("WignerSymbols")
Pkg.add("Memoize")
```

Then clone the repository:

```sh
git clone https://github.com/Seth-Kurankyi/su2bf-TNAlgo.git
```

### Usage
To use the functions provided in this repository, you can include the files in your Julia script or directly in the REPL. For example:

```julia
include("path_to_your_repository/coherent_vertex.jl")
include("path_to_your_repository/partial_cohn_vertex.jl")

result = coherent_vertex1([j12, j13, j14, j15, j23, j24, j25, j34, j35, j45], [nn1, nn2, nn3, nn4, nn5])
println(result)
```

### Contributions
Contributions are welcome! If you find any issues or have suggestions for improvements, please open an issue or submit a pull request.

### License
This project is licensed under the MIT License.

### References
If you use this code in your research, please cite the relevant papers and this repository.



