"""
Orbit determination functions.
"""

export gibbs, singleangle

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
