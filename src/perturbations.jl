"""
Functions for orbital perturbations.
"""

export geta, geti

"""
Given the mean perturbation in longitude of the ascending node, eccentricity, and inclination, estimate the semimajor axis.

# Arguments
- `barΩ`: mean perturbation in longitude of the ascending node, in degrees
- `e`: eccentricity
- `i`: inclination, in degrees

# Returns
- `a`: semimajor axis, in km
"""
function geta(barΩ::Float64, e::Float64, i::Float64)::Float64
    # Convert to radians
    barΩ = barΩ * (π / 180)
    i = i * (π / 180)

    # Compute some intermedial constants to keep equation clean
    k = (-3 / 2) * √(Earth.μ) * (Earth.J2 / (1 - e^2)^2) * (Earth.r)^2 * cos(i)

    a72 = k / barΩ

    return a72^(2 / 7)

end


"""
Given the mean perturbation in argument of perigee, eccentricity, and semimajor axis, estimate the inclination.

# Arguments
- `barω`: mean perturbation in argument of perigee, in degrees
- `e`: eccentricity
- `a`: semimajor axis, in km

# Returns
- `i`: inclination, in degrees
"""
function geti(barω::Float64, e::Float64, a::Float64)::Float64
    # Convert to radians
    barω = barω * (π / 180)

    # Compute some intermedial constants to keep equation clean
    n = √(Earth.μ / a^3)
    k = ((3 * n * Earth.J2) / (4 * (1 - e^2)^2)) * (Earth.r / a)^2

    sini = √((4 - (barω / k)) / 5)

    return asin(sini)

end