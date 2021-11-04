# DynamicFoam.jl

Generalizes [DynamicFoam](https://github.com/weigert/DynamicFoam) to 3D.

## Settings

Comparing to [DynamicFoam](https://github.com/weigert/DynamicFoam):
* The criticality used is equivalent to an infinite criticality and cannot be changed.
* If "Edge scale" is on, the dynamics is equivalent to [DynamicFoam](https://github.com/weigert/DynamicFoam): it scales the edges. If it is off, it scales delaunay simplices.
* The "Height" setting changes the height of the rendering frame but does not influence the dynamics. The points are are samples in a fixed square or cube that cannot be changed.

## Reporting issues

As the software use randomness. Running it twice give different results.
To make things reproducible, you can fix the seed as follows where 1234 can be replaced by any number:
```julia
julia> using DynamicFoam

julia> using Random

julia> Random.seed!(1234)

julia> foam(K, N)
```
If the program is run again with the same seed (1234 here), it will produce exactly the same simulations.

## Installation

Download this repository. Let `path/to/DynamicFoam` be the directory where this package was decompressed.

[Install Julia](https://julialang.org/downloads/) and launch it. You should see the following prompt:
```julia
julia>
```
Type `]` so that the prompts becomes
```julia
(@v1.6) pkg>
```
Typing backspace allows to go back to the prompt
```julia
julia>
```
Activate the environment of DynamicFoam.jl as follows:
```julia
(@v1.6) pkg> activate path/to/DynamicFoam
```
Now install its dependencies as follows:
```julia
(@v1.6) pkg> instantiate
```
The installation should take a while.
Once this is done, to launch the software with `K` points and `N` dimensions, do
```julia
julia> using DynamicFoam

julia> foam(K, N)
```
Note that `N` can only be 2 or 3.
For instance, for 10 points in 3D, do:
```julia
julia> foam(10, 3)
```

### Troubleshooting

On Arch Linux, linux should be launched as follows:
```julia
$ LD_PRELOAD=/usr/lib64/libstdc++.so.6 julia
```
See [this issue](https://github.com/JuliaGL/GLFW.jl/issues/198) for more details.
