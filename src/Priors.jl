using Distributions

"""
    mutable struct Priors
    
A mutable struct that encapsulates the prior distributions for various parameters in a Markov Chain Monte Carlo (MCMC) simulation.

# Fields   
- `logm1_dist`: Prior distribution for primary mass (in log10(Msun)).
- `logm2_dist`: Prior distribution for secondary mass (in log10(Msun)).
- `logP_dist`: Prior distribution for orbital period (in log10(days)).
- `vkick_dist`: Prior distribution for the kick velocity (in km/s).
- `frac_dist`: Prior distribution for the black hole formation fraction.
- `e_dist`: Prior distribution for the initial eccentricity.
- `Venv_α_100kms_dist`: Prior distribution for the environmnet RA velocity (in 100 km/s).
- `Venv_δ_100kms_dist`: Prior distribution for the environment Dec velocity (in 100 km/s).
- `Venv_r_100kms_dist`: Prior distribution for the environment radial velocity (in 100 km/s).
- `pm_α_dist`: Prior distribution for the proper motion in RA (in mas/yr).
- `pm_δ_dist`: Prior distribution for the proper motion in Dec (in mas/yr).
- `parallax_dist`: Prior distribution for the parallax (in mas).

Priors contains the prior distribution of each of the desired parameters. Many are not required, and some are 
degenerate with others. All default to missing. 

"""
@kwdef mutable struct Priors
    logm1_dist::Union{ContinuousUnivariateDistribution,Missing} = missing
    logm2_dist::Union{ContinuousUnivariateDistribution,Missing} = missing
    logP_dist::Union{ContinuousUnivariateDistribution,Missing} = missing 
    vkick_dist::Union{ContinuousUnivariateDistribution,Missing} = missing
    frac_dist::Union{ContinuousUnivariateDistribution,Missing} = missing 
    e_dist::Union{ContinuousUnivariateDistribution,Missing} = missing       
    Venv_α_100kms_dist::Union{ContinuousUnivariateDistribution,Missing} = missing
    Venv_δ_100kms_dist::Union{ContinuousUnivariateDistribution,Missing} = missing
    Venv_r_100kms_dist::Union{ContinuousUnivariateDistribution,Missing} = missing
    pm_α_dist::Union{ContinuousUnivariateDistribution,Missing} = missing  
    pm_δ_dist::Union{ContinuousUnivariateDistribution,Missing} = missing 
    parallax_dist::Union{ContinuousUnivariateDistribution,Missing} = missing 
end

"""
    Priors(keywords...)

Constructs a `Priors` object and its associated string representation, from a list of specified priors. 

# Arguments
- `keywords...`: prior distributions for various parameters.

# Returns
- A `Priors` object and its string representation. 

# Example
```
SideKicks.@Priors(
     logm1_dist = Uniform(0.1,3), 
     logm2_dist = Uniform(0.1,3), 
)
```
"""
macro Priors(keywords...)
    ex = Expr(:call)
    for keyword in keywords
        keyword.head = :kw
    end
    ex.args = [:(SideKicks.Priors), [keyword for keyword in keywords]...]
    return :($ex, $(string(ex)))
end
