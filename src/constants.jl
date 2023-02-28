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
end

const Earth = Planet(3.98600442e5, 6378.137)
