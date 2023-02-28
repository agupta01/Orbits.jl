"""
Orbit propagation functions. 
"""

export θ2E, E2θ, TOF, propagate, tof2RV, lambert

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
Compute lagrangian coefficients given r, r0, Δθ, and p.

# Arguments
- `r`: radius at t₁, in km
- `r0`: radius at t₀, in km
- `Δθ`: change in true anomaly, in radians
- `p`: parameter

# Returns
- `f`: lagrangian coefficient
- `g`: lagrangian coefficient
- `fdot`: lagrangian coefficient
- `gdot`: lagrangian coefficient
"""
function lagrangian(r::Float64, r0::Float64, Δθ::Float64, p::Float64)::Tuple{Float64,Float64,Float64,Float64}
    # Compute the lagrangian coefficients
    f = 1 - (r / p) * (1 - cos(Δθ))
    g = (r * r0 * sin(Δθ)) / (√(p * Earth.μ))
    fdot = √(Earth.μ / p) * (((1 - cos(Δθ)) / p) - (1 / r) - (1 / r0)) * tan(Δθ / 2)
    gdot = 1 - (r0 / p) * (1 - cos(Δθ))

    return f, g, fdot, gdot
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
    f, g, fdot, gdot = lagrangian(r, r0, Δθ, p)

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

"""
Solve Lambert's problem using the p-iteration method.

# Arguments
- `R1`: ECI position vector, in km at t₀
- `R2`: ECI position vector, in km at t₁
- `Δt`: time of flight, in minutes
- `longway`: boolean (Default false). If true, use the long way around the orbit. If false, use the short way.

# Returns
- `V1`: ECI velocity vector, in km/s at t₀
- `V2`: ECI velocity vector, in km/s at t₁
- `p`: Converged value of p, in km
- `e`: eccentricity
- `a`: semi-major axis, in km
"""
function lambert(R1::Vector{Float64}, R2::Vector{Float64}, Δt::Float64, longway::Bool=false)::Tuple{Vector{Float64},Vector{Float64},Float64,Float64,Float64}
    # Compute norms
    r1 = norm(R1)
    r2 = norm(R2)

    # Convert Δt to seconds
    Δt = Δt * 60

    # Compute Δθ
    Δθ = acos((R1 ⋅ R2) / (r1 * r2))
    if longway
        Δθ = 2π - Δθ
    end

    # Compute auxiliary components
    k = r1 * r2 * (1 - cos(Δθ))
    l = r1 + r2
    m = r1 * r2 * (1 + cos(Δθ))

    # calculate p's range
    pmin = k / (l + √(2m))
    pmax = k / (l - √(2m))

    # compute two initial guesses
    p1 = 0.7 * pmin + 0.3 * pmax
    p2 = 0.3 * pmin + 0.7 * pmax

    # compute the TOF of the second guess
    geta(p) = (m * k * p) / ((2m - l^2) * p^2 + (2 * k * l * p - k^2))
    function getTOF(p)
        f, g, fdot, gdot = lagrangian(r2, r1, Δθ, p)
        a = geta(p)
        cosΔE = 1 - (r1 / a) * (1 - f)
        sinΔE = -(r1 * r2 * fdot) / √(Earth.μ * a)
        ΔE = atan(sinΔE, cosΔE)
        if ΔE < 0
            ΔE += 2π
        end
        TOF = g + √(a^3 / Earth.μ) * (ΔE - sin(ΔE)) # in seconds
        return TOF
    end
    τ2 = getTOF(p1) - Δt
    τ = getTOF(p2) - Δt
    p = p2 - τ * ((p2 - p1) / (τ - τ2))
    while true
        # Compute the TOF of the new guess
        τ2 = τ
        τ = getTOF(p) - Δt
        println("p: $p, τ: $τ")
        (abs(τ) > 1e-4) || break # if τ is small enough, we're done

        # otherwise, compute the next guess
        p1 = p2
        p2 = p
        p = p2 - τ * ((p2 - p1) / (τ - τ2))
    end

    # Compute the Lagrangian coefficients to get V1 and V2
    f, g, fdot, gdot = lagrangian(r2, r1, Δθ, p)
    V1 = (1 / g) * (R2 - f * R1)
    V2 = fdot * R1 + (gdot / g) * (R2 - f * R1)

    # Compute the rest of the parameters
    a = geta(p)
    e = √(1 - p / a)

    return V1, V2, p, e, a
end
