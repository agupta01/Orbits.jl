# Orbits

A simple library for calculating various orbital mechanics parameters.

## Installation

```julia
julia> Pkg.add(PackageSpec(url="https://www.github.com/agupta01/Orbits.jl"))
```

## Usage

In this example, we'll get the apoapsis and periapsis for an earth orbiting satellite given current altitude, velocity, and flight path angle.
```julia
julia> using Orbits
julia> r = 1769.7560 + Earth.r # km
julia> v = 6.8164 # km/s
julia> γ = 4.9665 # degrees
julia> getapses(r, v, γ, "deg")
(6982.122754552764, 8533.768432997193)
```
