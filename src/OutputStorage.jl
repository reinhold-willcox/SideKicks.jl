using HDF5

"""
# OutputStorage.jl

This file provides functionality for saving and extracting MCMC results to and from HDF5 files. It includes:
- `SaveResults`: Saves MCMC results, metadata, and related information to an HDF5 file.
- `ExtractResults`: Extracts MCMC results, metadata, and related information from an HDF5 file.

These functions are designed to facilitate the storage and retrieval of simulation results in a structured and efficient manner.
"""

"""
    SaveResults(fname::String, kick_mcmc::KickMCMC)

Saves the results, metadata, and related information from a `KickMCMC` object to an HDF5 file.

# Arguments
- `fname::String`: The name of the HDF5 file to save the results to.
- `kick_mcmc::KickMCMC`: The `KickMCMC` object containing the results and metadata to save.

# Details
- Results are stored in a group named `"results"`.
- Observations and priors are stored as strings in a group named `"strings"`.
- Metadata, including NUTS sampler settings and summary statistics, is stored in a group named `"metadata"`.

# Example
```julia
SaveResults("results.h5", kick_mcmc)
```
""" 
function SaveResults(fname::String, kick_mcmc::KickMCMC) 
    h5open(fname, "w") do fid # Define results group 
        results = create_group(fid, "results") 
        dict_keys = keys(kick_mcmc.results) 
        results["results_keys"] = String.(dict_keys) 
        for key in dict_keys 
            results[String(key)] = kick_mcmc.results[key] 
        end

        # Define strings group
        strings = create_group(fid, "strings")
        strings["observations"] = kick_mcmc.observations_string
        strings["priors"] = kick_mcmc.priors_string

        # Define metadata group
        meta = create_group(fid, "metadata")
        meta["nuts_warmup_count"] = kick_mcmc.nuts_warmup_count
        meta["nuts_acceptance_rate"] = kick_mcmc.nuts_acceptance_rate
        meta["nsamples"] = kick_mcmc.nsamples

        # Add summary statistics in a subgroup of metadata
        stats = create_group(meta, "stats")
        ess = create_group(stats, "ess")
        dict_keys = keys(kick_mcmc.ess)
        ess["ess_keys"] = String.(dict_keys)
        for key in dict_keys
            ess[String(key)] = kick_mcmc.ess[key]
        end
        rhat = create_group(stats, "rhat")
        dict_keys = keys(kick_mcmc.rhat)
        rhat["rhat_keys"] = String.(dict_keys)
        for key in dict_keys
            rhat[String(key)] = kick_mcmc.rhat[key]
        end
    end
end

""" 
    ExtractResults(fname::String; transpose_results::Bool=false) -> Tuple

Extracts results, observations, priors, and metadata from an HDF5 file.

# Arguments
- `fname::String`: The name of the HDF5 file to extract results from.
- `transpose_results::Bool=false`: Whether to transpose 2D result arrays (default: false).

# Returns
- A tuple containing:

    1. `results::Dict{Symbol, Any}`: A dictionary of extracted results.
    2. `observations::Observations`: The observations object reconstructed from the file.
    3. `priors::Priors`: The priors object reconstructed from the file.
    4. `metadata::HDF5.Group`: The metadata group containing additional information.

# Details
- Results are extracted from the `"results"` group.
- Observations and priors are reconstructed from their string representations in the `"strings"` group.
- Metadata, including NUTS sampler settings and summary statistics, is extracted from the `"metadata"` group.

# Example
```julia
results, observations, priors, metadata = ExtractResults("results.h5")
```

# Throws
- `ArgumentError`: If the specified file does not exist. 
""" 
function ExtractResults(fname::String; transpose_results::Bool=false) 
    if ~isfile(fname) 
        throw(ArgumentError("File not found")) 
    end

    fid = h5open(fname, "r") 
    extracted_results = fid["results"] 
    results = Dict()

    # Extract results
    for key ∈ keys(extracted_results) 
        results[Symbol(key)] = extracted_results[key][] 
        if transpose_results && ndims(results[Symbol(key)]) == 2 
            results[Symbol(key)] = transpose(results[Symbol(key)]) 
        end 
    end

    # Extract strings
    strings = fid["strings"] 
    obs_string = strings["observations"][] 
    observations = eval(Meta.parse(obs_string)) 
    priors_string = strings["priors"][] 
    priors = eval(Meta.parse(priors_string))

    # Extract metadata
    metadata = fid["metadata"]
    return (results, observations, priors, metadata)
end

