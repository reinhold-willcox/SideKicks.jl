using StatsBase
using CairoMakie
using CornerPlotting

export create_corner_plot

"""
# CornerPlots.jl

This file provides functionality for creating corner plots, which are useful for visualizing the relationships between multiple parameters in MCMC (Markov Chain Monte Carlo) results. 

## Overview
- `PlottingProps`: A struct for storing plotting properties such as parameter names, units, ranges, and labels.
- `createPlottingProps`: A helper function to generate plotting properties from a matrix.
- `create_corner_plot`: The main function to generate corner plots for MCMC results.

These tools are designed to work with the `CornerPlotting` package and integrate with `CairoMakie` for high-quality visualizations.
"""

"""
    mutable struct PlottingProps

`PlottingProps` contains the properties, units, ranges, and names (LaTeX encouraged) for each parameter to plot.

# Fields
- `props::Vector{Symbol}`: The list of parameter names to plot.
- `units::Dict{Symbol, String}`: Units for each parameter.
- `ranges::Dict{Symbol, Tuple{Float64, Float64}}`: Ranges for each parameter.
- `names::Dict{Symbol, String}`: Display names (e.g., LaTeX labels) for each parameter.
"""
mutable struct PlottingProps
    props::Vector{Symbol}
    units::Dict{Symbol, String}
    ranges::Dict{Symbol, Tuple{Float64, Float64}}
    names::Dict{Symbol, String}
end

"""
    createPlottingProps(props_matrix::Vector{Vector{Any}}) -> Vector

Generates plotting properties from a matrix of parameter information.

# Arguments
- `props_matrix::Vector{Vector{Any}}`: A matrix where each row corresponds to a parameter and contains:
    - Name (Symbol)
    - Scaling factor (e.g., units)
    - Range (Tuple)
    - Label (String)

# Returns
A vector containing:
- `names`: A vector of parameter names.
- `scaling`: A dictionary mapping parameter names to scaling factors.
- `ranges`: A dictionary mapping parameter names to ranges.
- `labels`: A dictionary mapping parameter names to labels.
"""
function createPlottingProps(props_matrix::Vector{Vector{Any}})
    props_matrix = stack(props_matrix) # Convert to matrix
    names   = convert(Vector{Symbol}, props_matrix[1, :])
    scaling = Dict(zip(names, props_matrix[2, :]))
    ranges  = Dict(zip(names, props_matrix[3, :]))
    labels  = Dict(zip(names, props_matrix[4, :]))
    return [names, scaling, ranges, labels]
end

"""
    create_corner_plot(results, plotting_props; kwargs...) -> CornerPlot

Creates a corner plot for a selected subset of parameters from MCMC results.

# Arguments
- `results`: The extracted results (e.g., from an HDF5 object) from a previous MCMC run.
- `plotting_props::PlottingProps`: A `PlottingProps` object containing parameter properties and ranges.
- `dists_to_plot`: Optional distributions to overlay on the plots.
- `fig`: A `Figure` object for the plot (default: new figure).
- `supertitle`: A title for the entire plot (default: `nothing`).
- `fraction_1D`: The area fraction to include in the confidence interval bounds for 1D plots (default: `0.9`).
- `fractions_2D`: The area fractions to determine different colored regions in 2D plots (default: `[0.9]`).
- `show_CIs`: Whether to include confidence intervals (default: `true`).
- `nbins`: Number of bins for 1D histograms (default: `100`).
- `nbins_contour`: Number of bins for 2D contour plots (default: `30`).
- `supertitlefontsize`: Font size for the supertitle (default: `25`).
- `use_corner_plotting_theme`: Whether to use the default theme from `CornerPlotting` (default: `true`).

# Returns
- A `CornerPlot` object containing the generated corner plot.

# Notes
- If extra plotting properties are included that cannot be used, they will be ignored, and a warning will be printed.
- The function integrates with `CornerPlotting` for advanced plotting features.

# Example
```julia

corner_plot = create_corner_plot(results, plotting_props; supertitle="MCMC Results")
