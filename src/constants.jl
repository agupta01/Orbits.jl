"""
Constants (planet measurements, gravitational parameters, etc.)
"""

export Planet, Earth

"""
Models a planet.
"""
struct Planet
    μ::Float64 # km^3/s^2
    r::Float64 # km 
    J2::Float64 # J2 coefficient, measures the degree of oblateness of the planet
end

const Earth = Planet(3.98600442e5, 6378.137, 1.0826267e-3)


