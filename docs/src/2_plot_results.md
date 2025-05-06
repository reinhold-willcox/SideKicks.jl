```@meta
EditURL = "../../examples/plot_results_vfts243.jl"
```

# Plotting MCMC results from VFTS 243

In this example we produce corner plots from the results of VFTS 243 that were stored from
the previous example. We start by loading up

````@example 2_plot_results
using CairoMakie
using SideKicks
using CornerPlotting
using Distributions

results, observations, priors, metadata = SideKicks.ExtractResults("vfts243_results.hdf5")
````

As an initial check, we can plot the quantities that were used as observations. A basic
consistency check is to verify these are consistent with the MCMC samples.

````@example 2_plot_results
names = [:m1_f, :P_f, :e_f, :K1]
scaling = Dict( :m1_f => m_sun,
                :P_f  => day,
                :e_f  => 1,
                :K1   => km_per_s
               )
labels  = Dict( :m1_f => L"M_{1f}\;[M_{\odot}]",
                :P_f  => L"P_f\;[\mathrm{days}]",
                :e_f  => L"e_f",
                :K1   => L"K_1  \;[\mathrm{km s}^{-1}]",
               )
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

