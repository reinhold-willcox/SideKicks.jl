using Distributions

"""
    mutable struct Priors
    
TODO: Fix this
A mutable struct that encapsulates the prior distributions for various parameters in a Markov Chain Monte Carlo (MCMC) simulation.
# Fields   
- `logm1_dist`: Prior distribution for the logarithm of the first mass.
- `logm2_dist`: Prior distribution for the logarithm of the second mass.
- `logP_dist`: Prior distribution for the logarithm of the period.
- `vkick_dist`: Prior distribution for the kick velocity.
- `frac_dist`: Prior distribution for the fraction of stars kicked.
- `e_dist`: Prior distribution for the eccentricity.
- `Venv_α_100kms_dist`: Prior distribution for the alpha parameter of the velocity distribution at 100 km/s.
- `Venv_δ_100kms_dist`: Prior distribution for the delta parameter of the velocity distribution at 100 km/s.
- `Venv_r_100kms_dist`: Prior distribution for the radius parameter of the velocity distribution at 100 km/s.
- `rv_env_dist`: Prior distribution for the radial velocity of the environment.
- `pmra_dist`: Prior distribution for the proper motion in right ascension.
- `pmdec_dist`: Prior distribution for the proper motion in declination.
- `parallax_dist`: Prior distribution for the parallax.

Priors contains the prior distribution of each of the desired parameters

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
    rv_env_dist::Union{ContinuousUnivariateDistribution,Missing} = missing 
    pmra_dist::Union{ContinuousUnivariateDistribution,Missing} = missing  
    pmdec_dist::Union{ContinuousUnivariateDistribution,Missing} = missing 
    parallax_dist::Union{ContinuousUnivariateDistribution,Missing} = missing 
end

"""
    Priors(keywords...)

TODO: Fix this
Constructs a `Priors` object using keyword arguments.
# Arguments
- `keywords...`: Keyword arguments representing the prior distributions for various parameters.
# Returns
A `Priors` object initialized with the provided keyword arguments.
"""
macro Priors(keywords...)
    ex = Expr(:call)
    for keyword in keywords
        keyword.head = :kw
    end
    ex.args = [:(SideKicks.Priors), [keyword for keyword in keywords]...]
    return :($ex, $(string(ex)))
end
