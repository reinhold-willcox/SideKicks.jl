"""
# Observations.jl

This file defines the `Observations` struct and related functionality for handling observed data in the SideKicks.jl package.

## Overview
- `Observations`: A mutable struct for storing observed parameters, their values, errors, and units.
- `Observations(obs_matrix)`: A constructor function to create an `Observations` object from a matrix of observed data.
- `@Observations`: A macro for creating `Observations` objects programmatically.

These tools are designed to streamline the handling of observational data in astrophysical simulations.
"""

"""
    mutable struct Observations

`Observations` contains the symbols, values, errors, and units of each observed parameter.

# Fields
- `props::Vector{Symbol}`: A vector of parameter names (symbols) for the observed data.
- `vals::Vector{Float64}`: A vector of observed values for each parameter.
- `errs::Vector{Float64}`: A vector of errors associated with each observed value.
- `units::Vector{Float64}`: A vector of units for each observed parameter.
"""
@kwdef mutable struct Observations
    props::Vector{Symbol}
    vals::Vector{Float64}
    errs::Vector{Float64}
    units::Vector{Float64}
end

"""
    Observations(obs_matrix::Vector{Vector{Any}}) -> Observations

Constructs an `Observations` object from a matrix of observed data.

# Arguments
- `obs_matrix::Vector{Vector{Any}}`: A matrix where each row corresponds to a parameter and contains:
    - Name (Symbol)
    - Value (Float64)
    - Error (Float64)
    - Unit (Float64)

# Returns
An `Observations` object containing the observed parameters, their values, errors, and units.

# Example
```julia
obs_matrix = [
    [:param1, 1.0, 0.1, 1.0],
    [:param2, 2.0, 0.2, 1.0]
]
obs = Observations(obs_matrix)
""" 
function Observations(obs_matrix::Vector{Vector{Any}}) 
    obs_matrix = stack(obs_matrix) # Convert to a matrix 
    return Observations( 
        props = obs_matrix[1, :], 
        vals  = obs_matrix[2, :], 
        errs  = obs_matrix[3, :], 
        units = obs_matrix[4, :]) 
end

""" @Observations(observations)

A macro for programmatically creating Observations objects.

# Arguments
 - `observations`: A matrix or array-like structure containing observed data.

# Returns
- An expression that creates an Observations object.

# Example
```
obs = SideKicks.@Observations([
    [:P_f,  10.4031, 0.01,   day],
    [:e_f,  0.017,   0.012,  1],
])
```
""" 
macro Observations(observations) 
    ex = Expr(:call) 
    ex.args = [SideKicks.Observations, observations] 
    return :($ex, $(string(ex))) 
end
