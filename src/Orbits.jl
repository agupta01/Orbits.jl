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

"""
Determines a position vector given a single observation from a specified ground site.

# Arguments
- `ϕ`: latitude of the ground site, in degrees
- `λ`: longitude of the ground site, in degrees
- `ρ`: range from the ground site to the satellite, in km
- `β`: azimuth angle from the ground site to the satellite, in degrees
- `σ`: elevation angle from the ground site to the satellite, in degrees
- `angletype`: "deg" or "rad", depending on the units of the angles

# Returns
- `R`: position vector of the satellite in ECI coordinates, in km
"""
function singleangle(ϕ::Float64, λ::Float64, ρ::Float64, β::Float64, σ::Float64, angletype::String)::Vector{Float64}
    # Convert to radians if necessary
    if angletype == "deg"
        ϕ = ϕ * (π / 180)
        λ = λ * (π / 180)
        β = β * (π / 180)
        σ = σ * (π / 180)
    end

    I = [1.0, 0.0, 0.0]
    J = [0.0, 1.0, 0.0]
    K = [0.0, 0.0, 1.0]

    # Compute ρ_SEZ (range vector in SEZ coordinates)
    ρ_SEZ = ρ * (-cos(β) * cos(σ) * I + sin(β) * cos(σ) * J + sin(σ) * K)

    # Compute the transformation matrix
    sλ = sin(λ)
    cλ = cos(λ)
    sϕ = sin(ϕ)
    cϕ = cos(ϕ)
    D = zeros(3, 3)
    D[1, 1] = sϕ * cλ
    D[1, 2] = -sλ
    D[1, 3] = cϕ * cλ
    D[2, 1] = sϕ * sλ
    D[2, 2] = cλ
    D[2, 3] = cϕ * sλ
    D[3, 1] = -cϕ
    D[3, 3] = sϕ

    # Compute the range vector in ECI coordinates
    ρ_ECI = D * ρ_SEZ

    # Get the site vector
    R_site = Earth.r * (cϕ * cλ * I + cϕ * sλ * J + sϕ * K)

    return R_site + ρ_ECI
end

"""
Convert true anomaly to eccentric anomaly.

# Arguments
- `e`: eccentricity
- `θ`: true anomaly, in degrees
- `angletype`: "deg" or "rad", depending on the units of the angles

# Returns
- `E`: eccentric anomaly, in degrees
"""
function θ2E(e::Float64, θ::Float64, angletype::String)::Float64
    # Convert to radians if necessary
    if angletype == "deg"
        θ = θ * (π / 180)
    end

    # Compute the eccentric anomaly
    E = 2 * atan(√((1 - e) / (1 + e)) * tan(θ / 2))

    # Convert to degrees if necessary
    if angletype == "deg"
        E = E * (180 / π)
    end

    return E
end

"""
Convert eccentric anomaly to true anomaly.

# Arguments
- `e`: eccentricity
- `E`: eccentric anomaly, in degrees
- `angletype`: "deg" or "rad", depending on the units of the angles

# Returns
- `θ`: true anomaly, in degrees
"""
function E2θ(e::Float64, E::Float64, angletype::String)::Float64
    # Convert to radians if necessary
    if angletype == "deg"
        E = E * (π / 180)
    end

    # Compute the true anomaly
    θ = 2 * atan(√((1 + e) / (1 - e)) * tan(E / 2))

    # Convert to degrees if necessary
    if angletype == "deg"
        θ = θ * (180 / π)
    end

    return θ
end

"""
Determine time of flight given a (position, velocity) state and a true anomaly change.

# Arguments
- `R`: position vector, in km at t₀
- `V`: velocity vector, in km/s at t₀
- `Δθ`: change in true anomaly, in degrees

# Returns
- `Δt`: time of flight, in seconds
"""
function TOF(R::Vector{Float64}, V::Vector{Float64}, Δθ::Float64)::Float64
    # Compute the orbital elements at t₀
    a, e, i, Ω, ω, θ = RV2coe(R, V, "deg")

    # Compute the eccentric anomalies at t₀ and t₀ + Δt
    E0 = θ2E(e, θ, "deg")
    E1 = θ2E(e, θ + Δθ, "deg")

    # Fix E0 or E1 if they are negative
    if E0 < 0
        E0 += 360
    end
    if E1 < 0
        E1 += 360
    end

    # Convert to radians
    E0 = E0 * (π / 180)
    E1 = E1 * (π / 180)

    # Compute the mean motion
    n = √(Earth.μ / a^3)

    # Use Kepler's equation to get TOF
    Δt = (1 / n) * (E1 - e * sin(E1) - E0 + e * sin(E0))

    # If Δt is negative, add 2π/n
    if Δt < 0
        Δt += 2π / n
    end

    return Δt
end


"""
Propagate an orbit using Lagrangian coefficients.

# Arguments
- `R`: position vector, in km at t₀
- `V`: velocity vector, in km/s at t₀
- `Δθ`: change in true anomaly, in degrees

# Returns
- `R1`: position vector, in km at t₀ + Δt
- `V1`: velocity vector, in km/s at t₀ + Δt
"""
function propagate(R::Vector{Float64}, V::Vector{Float64}, Δθ::Float64)::Tuple{Vector{Float64},Vector{Float64}}
    # Compute the necessary orbital elements at t₀
    h = R × V
    p = norm(h)^2 / Earth.μ
    ξ = norm(V)^2 / 2 - Earth.μ / norm(R)
    a = -Earth.μ / (2 * ξ)
    e = √(1 - p / a)

    # Convert Δθ to radians
    Δθ = Δθ * (π / 180)

    # Determine initial and final true anomaly, in radians
    θ0 = acos((p / norm(R) - 1) / e)
    θ1 = θ0 + Δθ

    # Compute future r
    r0 = norm(R)
    r = p / (1 + e * cos(θ1))

    # Compute the Lagrangian coefficients
    f = 1 - (r / p) * (1 - cos(Δθ))
    g = (r * r0 * sin(Δθ)) / (√(p * Earth.μ))
    fdot = √(Earth.μ / p) * (((1 - cos(Δθ)) / p) - (1 / r) - (1 / r0)) * tan(Δθ / 2)
    gdot = 1 - (r0 / p) * (1 - cos(Δθ))

    # Calculate the new position and velocity
    R1 = f * R + g * V
    V1 = fdot * R + gdot * V

    return R1, V1
end


"""
Determine a satellite's new position and velocity given existing state and time of flight.

# Arguments
- `a`: semi-major axis, in km
- `e`: eccentricity
- `θ`: true anomaly, in degrees
- `Δt`: time of flight, in minutes

# Returns
- `r1`: propagated position radius, in km
- `v1`: propagated velocity magnitude, in km/s
- `γ1`: propagated flight path angle, in degrees
- `θ1`: propagated true anomaly, in degrees
"""
function tof2RV(a::Float64, e::Float64, θ::Float64, Δt::Float64)::Tuple{Float64,Float64,Float64,Float64}
    # Convert Δt to seconds and θ to radians
    Δt = Δt * 60
    θ = θ * (π / 180)

    # Compute the mean motion
    n = √(Earth.μ / a^3)

    # Compute the eccentric anomaly
    E = 2 * atan(√((1 - e) / (1 + e)) * tan(θ / 2))

    # Compute the mean anomaly
    M = E - e * sin(E)

    # Compute the propagated mean anomaly at t₁
    M1 = M + n * Δt

    # Use Newton's method to find a suitable E1 that satisfies Kepler's equation
    E1 = M1 + e * sin(M1) + (e^2 / 2) * sin(2 * M1)
    errfunc(E) = E - e * sin(E) - M1
    ∇errfunc(E) = 1 - e * cos(E)
    err = errfunc(E1)
    while abs(err) > 1e-6
        E1 = E1 - err / ∇errfunc(E1)
        err = errfunc(E1)
    end

    # Compute the propagated true anomaly at t₁
    θ1 = 2 * atan(√((1 + e) / (1 - e)) * tan(E1 / 2))

    # Compute the propagated radius at t₁
    r1 = a * (1 - e * cos(E1))

    # Compute the propagated velocity at t₁
    v1 = √(2 * Earth.μ / r1 - Earth.μ / a)

    # Compute the propagated flight path angle at t₁
    γ1 = atan((v1^2 / Earth.μ) * (1 + e * cos(θ1)) / (1 - e^2))

    # Convert to degrees
    θ1 = θ1 * (180 / π)
    γ1 = γ1 * (180 / π)

    return r1, v1, γ1, θ1
end

export Earth, getapses, RV2koe, koe2rv, RV2coe, coe2RV, makeR, gibbs, singleangle, θ2E, E2θ, TOF, propagate, tof2RV

end # module
