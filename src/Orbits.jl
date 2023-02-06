"""
A note on notation:

Vectors/Matrices will always be denoted with uppercase letters, e.g. R, V, H, etc.
Scalars will always be denoted with lowercase letters, e.g. r, v, h, etc.
"""
module Orbits

using LinearAlgebra

# Constants
struct Planet
    μ::Float64 # km^3/s^2
    r::Float64 # km 
end
const Earth = Planet(3.98600442e5, 6378.137)

# Functions
"""
Gives the apses of an earth orbit given the velocity, radius, and flight path angle.
Operates completely in the orbital plane (PQW/perifocal frame).

# Arguments
- `v`: velocity scalar, in km/s
- `r`: radius scalar, in km from center of orbiting body
- `γ`: flight path angle, in radians or degrees
- `angletype`: "deg" or "rad", depending on the units of the angles

# Returns
- `rp`: radius of periapsis, in km
- `ra`: radius of apoapsis, in km
"""
function getapses(v::Float64, r::Float64, γ::Float64, angletype::String)
    # Convert to radians if necessary
    if angletype == "deg"
        γ = γ * (π / 180)
    end

    μ = Earth.μ
    V = [v * sin(γ), v * cos(γ), 0]
    R = [r, 0, 0]
    E = (1 / μ) * ((v^2 - (μ / r)) * R - (R ⋅ V) * V)
    H = R × V
    e = norm(E)
    h = norm(H)
    p = h^2 / μ
    rp = p / (1 + e)
    ra = p / (1 - e)
    return rp, ra
end

"""
Computes the classical (Keplerian) orbit
elements given the Cartesian position and velocity vectors, the
gravitational parameter of the attracting body, and the desired angletype.
Operates completely in the orbital plane (PQW/perifocal frame).

# Arguments
- `R`: perifocal (PQW) position vector, in km
- `V`: perifocal (PQW) velocity vector, in km/s
- `angletype`: "deg" or "rad", depending on the desired units of the angles

# Returns
- `h`: specific angular momentum, in km^2/s
- `ξ`: specific energy, in km^2/s^2
- `a`: semimajor axis, in km
- `p`: parameter, in km
- `e`: eccentricity
- `T`: orbital period, in s
- `rp`: radius of periapsis, in km
- `ra`: radius of apoapsis, in km
- `γ`: flight path angle, in degrees or radians
- `θ`: true anomaly, in degrees or radians
"""
function RV2koe(R::Vector{Float64}, V::Vector{Float64}, angletype::String)
    μ = Earth.μ

    r = norm(R)
    v = norm(V)
    H = R × V
    h = norm(H)
    ξ = (v^2 / 2) - (μ / r)
    a = -μ / (2 * ξ)
    p = h^2 / μ
    e = √(1 - (p / a))
    rp = a * (1 - e)
    γ = atan(V[1] / V[2])
    θ = acos(((p / r) - 1) / e)

    if e >= 1
        # Hyperbolic or parabolic orbit
        ra = Inf
        T = nothing
    else
        # Elliptical orbit
        ra = a * (1 + e)
        T = (2π * √(a^3 / μ)) * (a^(3 / 2))
    end

    # Convert to degrees if desired
    if angletype == "deg"
        θ = θ * (180 / π)
        γ = γ * (180 / π)
    end

    return h, ξ, a, p, e, T, rp, ra, γ, θ
end

"""
Computes the position and velocity scalars and flight path angle
given the semimajor axis, eccentricity, and true anomaly (θ).

# Arguments
- `a`: semimajor axis, in km
- `e`: eccentricity
- `θ`: mean anomaly, in degrees or radians
- `angletype`: "deg" or "rad", depending on the units of the angles

# Returns
- `r`: radius, in km
- `v`: velocity, in km/s
- `γ`: flight path angle, in degrees or radians
"""
function koe2rv(a::Float64, e::Float64, θ::Float64, angletype::String)
    μ = Earth.μ

    # Convert to radians if necessary
    if angletype == "deg"
        θ = θ * (π / 180)
    end

    # Compute the radius
    p = a * (1 - e^2)
    r = p / (1 + e * cos(θ))

    # Compute the velocity
    ξ = -μ / (2 * a)
    v = √(2 * (ξ + (μ / r)))

    # Compute the flight path angle
    γ = atan(e * sin(θ) / (1 + e * cos(θ)))

    # Convert to degrees if necessary
    if angletype == "deg"
        γ = γ * (180 / π)
    end

    return r, v, γ
end

"""
Compute all 3D orbital elements given the cartesian position and velocity vectors.

# Arguments
- `R`: cartesian (ECI) position vector, in km
- `V`: cartesian (ECI) velocity vector, in km/s
- `angletype`: "deg" or "rad", depending on the desired units of the angles

# Returns
- `a`: semimajor axis, in km
- `e`: eccentricity
- `i`: inclination, in degrees or radians
- `Ω`: longitude of the ascending node, in degrees or radians
- `ω`: argument of periapsis, in degrees or radians
- `θ`: true anomaly, in degrees or radians
"""
function RV2coe(R::Vector{Float64}, V::Vector{Float64}, angletype::String)
    μ = Earth.μ

    # ECI unit vectors
    I = [1.0, 0.0, 0.0]
    J = [0.0, 1.0, 0.0]
    K = [0.0, 0.0, 1.0]

    # R and V scalars
    r = norm(R)
    v = norm(V)

    # Semimajor axis
    ξ = (v^2 / 2) - (μ / r)
    a = -μ / (2 * ξ)

    # Eccentricity
    E = (1 / μ) * ((v^2 - (μ / r)) * R - (R ⋅ V) * V)
    e = norm(E)

    # Inclination
    H = R × V
    h = norm(H)
    i = acos(K ⋅ H / h)

    # Longitude of the ascending node
    N = K × H
    n = norm(N)
    cosΩ = (I ⋅ N) / n
    sinΩ = (J ⋅ N) / n
    Ω = atan(sinΩ, cosΩ)

    # Argument of periapsis
    if E[3] >= 0
        ω = acos(N ⋅ E / (n * e))
    else
        ω = 2π - acos(N ⋅ E / (n * e))
    end

    # True anomaly
    if (R ⋅ V) >= 0
        θ = acos(E ⋅ R / (e * r))
    else
        θ = 2π - acos(E ⋅ R / (e * r))
    end

    if isapprox(i, 0)
        if isapprox(e, 0)
            # Circular equatorial orbit
            θ = Ω + ω + θ # replace with true longitude at epoch
            Ω = nothing
            ω = nothing
        else
            # Elliptical equatorial orbit
            ω = Ω + ω # replace with longitude of periapsis
            Ω = nothing
        end
    end

    # Convert to degrees if necessary
    if angletype == "deg"
        i = i * (180 / π)
        Ω = Ω * (180 / π)
        ω = ω * (180 / π)
        θ = θ * (180 / π)
    end

    return a, e, i, Ω, ω, θ
end

"""
Compute the position and velocity vectors given the orbital elements.

# Arguments
- `a`: semimajor axis, in km
- `e`: eccentricity
- `i`: inclination, in degrees or radians
- `Ω`: longitude of the ascending node, in degrees or radians
- `ω`: argument of periapsis, in degrees or radians
- `θ`: true anomaly, in degrees or radians
- `angletype`: "deg" or "rad", depending on the units of the angles

# Returns
- `R`: cartesian (ECI) position vector, in km
- `V`: cartesian (ECI) velocity vector, in km/s
"""
function coe2RV(a::Float64, e::Float64, i::Float64, Ω::Float64, ω::Float64, θ::Float64, angletype::String)
    μ = Earth.μ

    # Convert to radians if necessary
    if angletype == "deg"
        i = i * (π / 180)
        Ω = Ω * (π / 180)
        ω = ω * (π / 180)
        θ = θ * (π / 180)
    end

    # Compute preliminaries
    p = a * (1 - e^2)
    r = p / (1 + e * cos(θ))
    h = √(p * μ)

    sθ = sin(θ)
    cθ = cos(θ)

    # Compute the position and velocity vectors in PQW
    R = r * [cθ, sθ, 0.0]
    V = (μ / h) * [-sθ, e + cθ, 0.0]

    # Print the PQW vectors
    println("PQW vectors:")
    println("R = $(R)")
    println("V = $(V)")

    # Rotate the vectors into the ECI frame
    R_transform = makeR(i, Ω, ω, "rad")
    R = R_transform * R
    V = R_transform * V

    return R, V
end

"""
Create a transformation matrix to convert PQW coordinates to ECI coordinates.

# Arguments
- `i`: inclination, in degrees or radians
- `Ω`: longitude of the ascending node, in degrees or radians
- `ω`: argument of periapsis, in degrees or radians
- `angletype`: "deg" or "rad", depending on the units of the angles

# Returns
- `R`: transformation matrix
"""
function makeR(i::Float64, Ω::Float64, ω::Float64, angletype::String)
    # Convert to radians if necessary
    if angletype == "deg"
        i = i * (π / 180)
        Ω = Ω * (π / 180)
        ω = ω * (π / 180)
    end

    sω = sin(ω)
    cω = cos(ω)
    sΩ = sin(Ω)
    cΩ = cos(Ω)
    si = sin(i)
    ci = cos(i)

    R = zeros(3, 3)

    R[1, 1] = (cΩ * cω) - (sΩ * sω * ci)
    R[1, 2] = (-cΩ * sω) - (sΩ * cω * ci)
    R[1, 3] = sΩ * si
    R[2, 1] = (sΩ * cω) + (cΩ * sω * ci)
    R[2, 2] = (-sΩ * sω) + (cΩ * cω * ci)
    R[2, 3] = -cΩ * si
    R[3, 1] = sω * si
    R[3, 2] = cω * si
    R[3, 3] = ci

    return R
end


"""
Uses Gibbs' method of orbit determination to compute the satellite's
intertial velocity vector V2 based on 3 cartesian (ECI) position 
vectors: R1, R2, R3.

# Arguments
- `R1`: cartesian (ECI) position vector at time t1, in km
- `R2`: cartesian (ECI) position vector at time t2, in km
- `R3`: cartesian (ECI) position vector at time t3, in km
- `getcoe`: if true, return the orbital elements instead of the V2 velocity vector

# Returns
if `getcoe` is false:
- `V2`: cartesian (ECI) velocity vector at time t2, in km/s
if `getcoe` is true:
- `a`: semimajor axis, in km
- `e`: eccentricity
- `i`: inclination, in degrees
- `Ω`: longitude of the ascending node, in degrees
- `ω`: argument of periapsis, in degrees
- `θ`: true anomaly, in degrees
"""
function gibbs(R1::Vector{Float64}, R2::Vector{Float64}, R3::Vector{Float64}, getcoe::Bool=false)::Union{Vector{Float64},Tuple{Float64,Float64,Float64,Float64,Float64,Float64}}
    μ = Earth.μ

    # Compute norms
    r1 = norm(R1)
    r2 = norm(R2)
    r3 = norm(R3)

    # Check if position vectors are coplanar
    ϵ = (R1 ⋅ (R2 × R3)) / (r1 * norm(R2 × R3))
    if abs(ϵ) >= 0.0349
        @error "Position vectors are not coplanar."
    end

    # Compute angles
    θ12 = acos(R1 ⋅ R2 / (r1 * r2))
    θ23 = acos(R2 ⋅ R3 / (r2 * r3))

    if isapprox(θ12, 0.0, atol=(π / 180)) || isapprox(θ23, 0.0, atol=(π / 180))
        @warn "Position vectors are too close together. Results may be inaccurate."
    end

    # Compute auxiliary vectors
    D = (R2 × R3) + (R3 × R1) + (R1 × R2)
    N = r1 * (R2 × R3) + r2 * (R3 × R1) + r3 * (R1 × R2)
    S = (r2 - r3) * R1 + (r3 - r1) * R2 + (r1 - r2) * R3

    # Compute the velocity vector
    V2 = (1 / r2) * √(μ / (norm(N) * norm(D))) * (D × R2) + √(μ / (norm(N) * norm(D))) * S

    if getcoe
        # Compute the orbital elements
        return RV2coe(R2, V2, "deg")
    else
        return V2
    end
end


export Earth, getapses, RV2koe, koe2rv, RV2coe, coe2RV, makeR, gibbs

end # module
