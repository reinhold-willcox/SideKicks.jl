using Turing
using Distributions

"""
# TuringModels.jl

This file defines custom probability distributions and Turing models for MCMC sampling of pre- and post-explosion properties of astrophysical systems.

## Overview
- `WrappedCauchy`: A custom distribution resembling the Cauchy distribution on the unit circle.
- `ShiftedVonMises`: A wrapper for the VonMises distribution to extend its domain to [0, 2π].
- `create_simplified_mcmc_model`: Creates a simplified Turing model for MCMC sampling, assuming pre-explosion circularity.
- `create_general_mcmc_model`: Creates a general Turing model for MCMC sampling, allowing for pre-explosion eccentricity.
"""

"""
    struct WrappedCauchy{T1<:Real, T2<:Real} <: ContinuousUnivariateDistribution

The `WrappedCauchy` distribution resembles the Cauchy distribution defined on the unit circle from 0 to 2π, with the endpoints wrapped back to each other.

# Fields
- `μ::T1`: The location parameter (mean angle).
- `σ::T2`: The scale parameter (controls the spread).
"""
struct WrappedCauchy{T1<:Real, T2<:Real} <: ContinuousUnivariateDistribution
    μ::T1
    σ::T2
end

"""
TODO: add docs
"""
Distributions.logpdf(d::WrappedCauchy, x::Real) = log(1 / (2 * π) * sinh(d.σ)) - log(cosh(d.σ) - cos(x - d.μ))
Distributions.pdf(d::WrappedCauchy, x::Real) = 1 / (2 * π) * sinh(d.σ) / (cosh(d.σ) - cos(x - d.μ))







"""
    struct ShiftedVonMises{T1<:Real, T2<:Real} <: ContinuousUnivariateDistribution

The `ShiftedVonMises` distribution is a wrapper for the VonMises distribution (as defined in Distributions.jl) to extend its domain to [0, 2π].

# Fields
- `μ::T1`: The location parameter (mean angle).
- `κ::T2`: The concentration parameter (controls the spread).
- `vonMisesDist::ContinuousUnivariateDistribution`: The underlying VonMises distribution.
"""
struct ShiftedVonMises{T1<:Real, T2<:Real, T3<:ContinuousUnivariateDistribution} <: ContinuousUnivariateDistribution
    μ::T1
    κ::T2
    vonMisesDist::T3
end

"""
    ShiftedVonMises(μ::Real, κ::Real) -> ShiftedVonMises

Creates a `ShiftedVonMises` distribution, which is a wrapper for the VonMises distribution to extend its domain to [0, 2π].

# Arguments
- `μ::Real`: The location parameter (mean angle) of the distribution.
- `κ::Real`: The concentration parameter (controls the spread) of the distribution.

# Returns
A `ShiftedVonMises` object with the specified parameters.

# Notes
- The `ShiftedVonMises` distribution ensures that the VonMises distribution is properly wrapped around the unit circle, making it suitable for angular data.
"""
function wrapped_vonmises(μ, κ)
    return ShiftedVonMises(μ, κ, VonMises(μ, κ))
end

"""
    Distributions.logpdf(d::ShiftedVonMises, x::Real) -> Real

Computes the log-probability density function (log-PDF) of the `ShiftedVonMises` distribution at a given value `x`.

# Arguments
- `d::ShiftedVonMises`: The `ShiftedVonMises` distribution object.
- `x::Real`: The value at which to evaluate the log-PDF.

# Returns
The log-probability density of the `ShiftedVonMises` distribution at `x`.

# Notes
- If `x` is outside the range `[μ - π, μ + π]`, it is shifted by an appropriate multiple of `2π` to bring it within this range before evaluating the log-PDF.
- This ensures the periodicity of the distribution on the unit circle.
"""
function Distributions.logpdf(d::ShiftedVonMises, x::Real)
    # If x is outside the range [μ-π, μ+π], we need to shift it by the correct
    # amount of 2π to fit it there
    if x > d.μ + π
        return logpdf(d.vonMisesDist, x - 2π)
    elseif x < d.μ - π
        return logpdf(d.vonMisesDist, x + 2π)
    else
        return logpdf(d.vonMisesDist, x)
    end
end

"""
    Distributions.pdf(d::ShiftedVonMises, x::Real) -> Real

Computes the probability density function (PDF) of the `ShiftedVonMises` distribution at a given value `x`.

# Arguments
- `d::ShiftedVonMises`: The `ShiftedVonMises` distribution object.
- `x::Real`: The value at which to evaluate the PDF.

# Returns
The probability density of the `ShiftedVonMises` distribution at `x`.

# Notes
- If `x` is outside the range `[μ - π, μ + π]`, it is shifted by an appropriate multiple of `2π` to bring it within this range before evaluating the PDF.
- This ensures the periodicity of the distribution on the unit circle.
"""
function Distributions.pdf(d::ShiftedVonMises, x::Real)
    # If x is outside the range [μ-π, μ+π], we need to shift it by the correct
    # amount of 2π to fit it there
    if x > d.μ + π
        return pdf(d.vonMisesDist, x - 2π)
    elseif x < d.μ - π
        return pdf(d.vonMisesDist, x + 2π)
    else
        return pdf(d.vonMisesDist, x)
    end
end

"""
    create_simplified_mcmc_model(; observations, priors, likelihood=:Cauchy, bhModel=arbitraryEjectaBH)

Creates a Turing model to perform a simplified MCMC sampling of the pre-explosion and kick properties of a system, assuming pre-explosion circularity.

# Arguments
- `observations::Observations`: Observational data, including parameters, values, errors, and units.
- `priors::Priors`: Prior distributions for the parameters of interest.
- `likelihood`: The likelihood distribution to use (`:Cauchy` or `:Normal`).
- `bhModel`: The model for black hole formation (default: `arbitraryEjectaBH`).

# Returns
- A tuple containing the Turing model and the list of return properties.
"""
function create_simplified_mcmc_model(;
    observations::Observations,
    priors::Priors,
    likelihood = :Cauchy,
    bhModel = arbitraryEjectaBH
    )

    logm1_dist = priors.logm1_dist
    logm2_dist = priors.logm2_dist
    logP_dist = priors.logP_dist
    vkick_dist = priors.vkick_dist
    frac_dist = priors.frac_dist

    valid_values = [:P_f, :P_circ, :e_f, :K1, :K2, :m1_f, :m2_f, :i_f, :Δv]
    for prop ∈ observations.props
        if prop ∉ valid_values
            throw(DomainError(observations.props, "Allowed observations are only [:P_f, :e_f, :K1, :K2, :m1_f, :m2_f, :i_f, :Δv]"))
        end
    end
    use_Pf = true
    if :P_circ ∈ observations.props
        if (:P_f ∈ observations.props) || (:e_f ∈ observations.props)
            throw(DomainError(observations.props, "Cannot have P_circ if either P_f or e_f are included"))
        else
            use_Pf = false
        end
    end

    if !(likelihood == :Cauchy || likelihood == :Normal)
        throw(DomainError(likelihood, "likelihood must be either :Cauchy or :Normal"))
    end

    @model function create_mcmc_model(props, obs_vals, obs_errs) 

        # set priors
        #Pre-explosion masses and orbital period
        logm1 ~ logm1_dist
        m1_i = 10^(logm1)*m_sun
        logm2 ~ logm2_dist
        m2_i = 10^(logm2)*m_sun
        logP ~ logP_dist
        P_i = 10^(logP)*day
        a_i = kepler_a_from_P(m1=m1_i, m2=m2_i, P=P_i)
        cosi ~ Uniform(0,1)
        i_f = acos(cosi)

        #Post-explosion masses
        frac ~ frac_dist
        m2_f = bhModel(m2_i, frac) # star 2 explodes, star 1 is kept fixed

        #Kick parameters
        vkick_100kms ~ vkick_dist  
        vkick = vkick_100kms *100*km_per_s
        cosθ ~ Uniform(-1,1)
        θ = acos(cosθ)
        xϕ ~ Normal(0,1)
        yϕ ~ Normal(0,1)
        normϕ = 1/sqrt(xϕ^2 + yϕ^2)
        cosϕ = xϕ*normϕ
        ϕ = acos(cosϕ)
        if yϕ < 0
            ϕ = 2π - ϕ
        end

        #m1 is assumed to remain constant
        m1_f = m1_i
        a_f, e_f = post_supernova_circular_orbit_a(m1_i=m1_i, m2_i=m2_i, a_i=a_i, m2_f=m2_f, vkick=vkick, θ=θ, ϕ=ϕ)
        P_f = kepler_P_from_a(m1=m1_f, m2=m2_f, a=a_f)
        P_circ = P_f * (1-e_f^2)^(3/2)
        if use_Pf 
            K1 = RV_semiamplitude_K1(m1=m1_f, m2=m2_f, P=P_f, e=e_f, i=i_f)
            K2 = RV_semiamplitude_K1(m1=m2_f, m2=m1_f, P=P_f, e=e_f, i=i_f)
        else
            K1 = RV_semiamplitude_K1(m1=m1_f, m2=m2_f, P=P_circ, e=0.0, i=i_f)
            K2 = RV_semiamplitude_K1(m1=m2_f, m2=m1_f, P=P_circ, e=0.0, i=i_f)
        end
        Δv = post_supernova_circular_orbit_Δv(m1_i=m1_i, m2_i=m2_i, a_i=a_i, m1_f=m1_i, m2_f=m2_f, vkick=vkick, θ=θ, ϕ=ϕ)

        use_cauchy = likelihood == :Cauchy
        for ii in eachindex(props)
            obs_symbol = props[ii]
            if obs_symbol == :P_f
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(P_f, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(P_f, obs_errs[ii])
            elseif obs_symbol == :P_circ
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(P_circ, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(P_circ, obs_errs[ii])
            elseif obs_symbol == :e_f
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(e_f, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(e_f, obs_errs[ii])
            elseif obs_symbol == :K1
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(K1, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(K1, obs_errs[ii])
            elseif obs_symbol == :K2
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(K2, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(K2, obs_errs[ii])
            elseif obs_symbol == :m1_f
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(m1_f, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(m1_f, obs_errs[ii])
            elseif obs_symbol == :m2_f
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(m2_f, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(m2_f, obs_errs[ii])
            elseif obs_symbol == :i_f
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(i_f, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(i_f, obs_errs[ii])
            elseif obs_symbol == :Δv
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(Δv, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(Δv, obs_errs[ii])
            end
        end

        # other params
        dm2 = m2_i - m2_f
        return     ( m1_i,    m2_i,    P_i,   a_i,     i_f,  
                     vkick,  θ,  ϕ,  m2_f,  a_f,   P_f,  P_circ,  e_f,  
                     K1,  K2,  frac,  dm2,  Δv) 
    end
    return_props = [:m1_i,   :m2_i,   :P_i,  :a_i,    :i_f, 
                    :vkick, :θ, :ϕ, :m2_f, :a_f,  :P_f, :P_circ, :e_f, 
                    :K1, :K2, :frac, :dm2, :Δv] 

    obs_vals_cgs = observations.vals .* observations.units
    obs_errs_cgs = observations.errs .* observations.units
    return [create_mcmc_model(observations.props, obs_vals_cgs, obs_errs_cgs), return_props]
end

"""
    create_general_mcmc_model(; observations, priors, likelihood=:Cauchy, bhModel=arbitraryEjectaBH)

Creates a Turing model to perform an MCMC sampling of the pre-explosion and kick properties of a system, allowing for pre-explosion eccentricity.

# Arguments
- `observations::Observations`: Observational data, including parameters, values, errors, and units.
- `priors::Priors`: Prior distributions for the parameters of interest.
- `likelihood`: The likelihood distribution to use (`:Cauchy` or `:Normal`).
- `bhModel`: The model for black hole formation (default: `arbitraryEjectaBH`).

# Returns
- A tuple containing the Turing model and the list of return properties.
"""
function create_general_mcmc_model(;
    observations::Observations,
    priors::Priors,
    likelihood = :Cauchy,
    bhModel = arbitraryEjectaBH
    )
    
    logm1_dist = priors.logm1_dist
    logm2_dist = priors.logm2_dist
    logP_dist = priors.logP_dist
    e_dist = priors.e_dist
    vkick_dist = priors.vkick_dist
    frac_dist = priors.frac_dist
    Venv_α_100kms_dist = priors.Venv_α_100kms_dist
    Venv_δ_100kms_dist = priors.Venv_δ_100kms_dist
    Venv_r_100kms_dist = priors.Venv_r_100kms_dist

    valid_values = [:P_f, :P_circ, :e_f, :K1, :K2, :m1_f, :m2_f, :Ω_f, :ω_f, :i_f, :vf_α, :vf_δ, :vf_r]
    for prop ∈ observations.props
        if prop ∉ valid_values
            throw(DomainError(observations.props, "Allowed observations are only [:P_f, :P_circ, :e_f, :K1, :K2, :m1_f, :m2_f, :Ω_f, :ω_f, :i_f, :vf_α, :vf_δ, :vf_r]"))
        end
    end
    use_Pf = true
    if :P_circ ∈ observations.props
        if (:P_f ∈ observations.props) || (:e_f ∈ observations.props)
            throw(DomainError(observations.props, "Cannot have P_circ if either P_f or e_f are included"))
        else
            use_Pf = false
        end
    end

    if !(likelihood == :Cauchy || likelihood == :Normal)
        throw(DomainError(likelihood, "likelihood must be either :Cauchy or :Normal"))
    end

    @model function create_mcmc_model(props, obs_vals, obs_errs) 
        # set priors
        #Pre-explosion masses and orbital period
        logm1 ~ logm1_dist
        m1_i = 10^(logm1)*m_sun
        logm2 ~ logm2_dist
        m2_i = 10^(logm2)*m_sun
        logP ~ logP_dist
        P_i = 10^(logP)*day
        a_i = kepler_a_from_P(m1=m1_i, m2=m2_i, P=P_i)
        e_i ~ e_dist

        #azimuthal angles are computed by sampling random points with a circularly symmetric distribution
        #we take all true anomalies to be equally likely. This is corrected by weighting later
        xν ~ Normal()
        yν ~ Normal()
        normν = 1/sqrt(xν^2+yν^2)
        cosν = xν*normν
        ν_i = acos(cosν)
        if yν < 0
            ν_i = 2π - ν_i
        end
        xΩ ~ Normal()
        yΩ ~ Normal()
        normΩ = 1/sqrt(xΩ^2+yΩ^2)
        cosΩ = xΩ*normΩ
        Ω_i = acos(cosΩ)
        if yΩ < 0
            Ω_i = 2π - Ω_i
        end
        xω ~ Normal()
        yω ~ Normal()
        normω = 1/sqrt(xω^2+yω^2)
        cosω = xω*normω
        ω_i = acos(cosω)
        if yω < 0
            ω_i = 2π - ω_i
        end
        cosi ~ Uniform(0,1)
        i_i = acos(cosi)

        #Post-explosion masses
        frac ~ frac_dist
        m2_f = bhModel(m2_i, frac) # star 2 explodes, star 1 is kept fixed

        #Kick parameters
        vkick_100kms ~ vkick_dist
        vkick = vkick_100kms*100*km_per_s
        cosθ ~ Uniform(-1,1)
        θ = acos(cosθ)
        xϕ ~ Normal(0,1)
        yϕ ~ Normal(0,1)
        normϕ = 1/sqrt(xϕ^2+yϕ^2)
        cosϕ = xϕ*normϕ
        ϕ = acos(cosϕ)
        if yϕ < 0
            ϕ = 2π - ϕ
        end

        #Initial systemic velocity parameters
        vi_α_100kms ~ Venv_α_100kms_dist
        vi_δ_100kms ~ Venv_δ_100kms_dist
        vi_r_100kms ~ Venv_r_100kms_dist
        vi_α = vi_α_100kms*100*km_per_s 
        vi_δ = vi_δ_100kms*100*km_per_s 
        vi_r = vi_r_100kms*100*km_per_s 

        #m1 is assumed to remain constant, no impact velocity - TODO: relax this later
        m1_f = m1_i
        a_f, e_f, Ω_f, ω_f, i_f, Δv_α, Δv_δ, Δv_r = 
            post_supernova_general_orbit_parameters(m1_i=m1_i, m2_i=m2_i, a_i=a_i, e_i=e_i, m2_f=m2_f, 
                vkick=vkick, θ=θ, ϕ=ϕ, ν_i=ν_i, Ω_i=Ω_i, ω_i=ω_i, i_i=i_i)
        P_f = kepler_P_from_a(m1=m1_f, m2=m2_f, a=a_f)
        P_circ = P_f * (1-e_f^2)^(3/2)
        if use_Pf 
            K1 = RV_semiamplitude_K1(m1=m1_f, m2=m2_f, P=P_f, e=e_f, i=i_f)
            K2 = RV_semiamplitude_K1(m1=m2_f, m2=m1_f, P=P_f, e=e_f, i=i_f)
        else
            K1 = RV_semiamplitude_K1(m1=m1_f, m2=m2_f, P=P_circ, e=0.0, i=i_f)
            K2 = RV_semiamplitude_K1(m1=m2_f, m2=m1_f, P=P_circ, e=0.0, i=i_f)
        end
        Δv = sqrt( Δv_δ^2 + Δv_α^2 + Δv_r^2)

        vf_α = vi_α + Δv_α
        vf_δ = vi_δ + Δv_δ
        vf_r = vi_r + Δv_r
        
        use_cauchy = likelihood == :Cauchy
        for ii in eachindex(props)
            obs_symbol = props[ii]
            if obs_symbol == :P_f
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(P_f, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(P_f, obs_errs[ii])
            elseif obs_symbol == :P_circ
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(P_circ, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(P_circ, obs_errs[ii])
            elseif obs_symbol == :e_f
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(e_f, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(e_f, obs_errs[ii])
            elseif obs_symbol == :K1
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(K1, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(K1, obs_errs[ii])
            elseif obs_symbol == :K2
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(K2, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(K2, obs_errs[ii])
            elseif obs_symbol == :m1_f
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(m1_f, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(m1_f, obs_errs[ii])
            elseif obs_symbol == :m2_f
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(m2_f, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(m2_f, obs_errs[ii])
            elseif obs_symbol == :i_f
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(i_f, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(i_f, obs_errs[ii])
            elseif obs_symbol == :Ω_f
                use_cauchy ?
                    obs_vals[ii] ~ WrappedCauchy(Ω_f, obs_errs[ii]) :
                    obs_vals[ii] ~ ShiftedVonMises(Ω_f, 1/obs_errs[ii]^2)
            elseif obs_symbol == :ω_f
                use_cauchy ?
                    obs_vals[ii] ~ WrappedCauchy(ω_f, obs_errs[ii]) :
                    obs_vals[ii] ~ ShiftedVonMises(ω_f, 1/obs_errs[ii]^2)
            elseif obs_symbol == :vf_α
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(vf_α, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(vf_α, obs_errs[ii])
            elseif obs_symbol == :vf_δ
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(vf_δ, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(vf_δ, obs_errs[ii])
            elseif obs_symbol == :vf_r
                use_cauchy ?
                    obs_vals[ii] ~ Cauchy(vf_r, obs_errs[ii]) :
                    obs_vals[ii] ~ Normal(vf_r, obs_errs[ii])
            end
        end

        dm2 = m2_i - m2_f
        sum_ωi_νi = (ω_i + ν_i) % 2π # should be between 0 and 2π

        return     ( m1_i,  m2_i,  P_i,  e_i,  a_i,  Ω_i,  ω_i,  i_i,  sum_ωi_νi,
                     vkick,  θ,  ϕ,  dm2,  frac,  ν_i, 
                     m1_f,  m2_f,  P_f,  P_circ, e_f,  a_f,  Ω_f,  ω_f,  i_f, 
                     K1,  K2,  vi_α,  vi_δ,  vi_r,  vf_α,  vf_δ,  vf_r, 
                     Δv_α,  Δv_δ,  Δv_r,  Δv )
    end
    return_props = [:m1_i, :m2_i, :P_i, :e_i, :a_i, :Ω_i, :ω_i, :i_i, :sum_ωi_νi,
                    :vkick, :θ, :ϕ, :dm2, :frac, :ν_i, 
                    :m1_f, :m2_f, :P_f, :P_circ, :e_f, :a_f, :Ω_f, :ω_f, :i_f, 
                    :K1, :K2, :vi_α, :vi_δ, :vi_r, :vf_α, :vf_δ, :vf_r, 
                    :Δv_α, :Δv_δ, :Δv_r, :Δv ]

    # Need to combine some of the observations to compare against the predicted output
    obs_vals_cgs = observations.vals .* observations.units
    obs_errs_cgs = observations.errs .* observations.units

    return [create_mcmc_model(observations.props, obs_vals_cgs, obs_errs_cgs), return_props]
    
end
