"""
A note on notation:

Vectors/Matrices will always be denoted with uppercase letters, e.g. R, V, H, etc.
Scalars will always be denoted with lowercase letters, e.g. r, v, h, etc.
"""
module Orbits

using LinearAlgebra

# Constants
include("./constants.jl")

# Functions
include("./parameters.jl") # ch. 2
include("./determination.jl") # ch. 3
include("./propagation.jl") # ch. 4
include("./perturbations.jl") # ch. 5

end # module
