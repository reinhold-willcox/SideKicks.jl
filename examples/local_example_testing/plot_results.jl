#=
# Plotting MCMC results 
=#

using CairoMakie
using SideKicks
using CornerPlotting

pwd = String(@__DIR__)*"/"
results, observations, priors, metadata = SideKicks.ExtractResults(pwd*"results.hdf5")

system_name = "Example 1: e_i=0"
system_fname = pwd*"example1"

##

plotting_props_obs_check = SideKicks.createPlottingProps([
   [:m1_f,    m_sun,     [10,20], L"M_1\;[M_{\odot}]"],
   [:m2_f,    m_sun,     [0, 20], L"M_2  \;[M_{\odot}]"],
   [:P_f,     day,       missing, L"P_f\;[\mathrm{days}]"],
   [:e_f,     1,         [0.1, 0.3], L"e_f"],
   [:K1,      km_per_s,  [0, 100], L"K_1  \;[\mathrm{km s}^{-1}]"],
   #[:v_N,     km_per_s, [130,170],      L"v_N  \;[\mathrm{km s}^{-1}]"],
   #[:v_E,     km_per_s, [380,430],      L"v_E  \;[\mathrm{km s}^{-1}]"],
   #[:v_r,     km_per_s, [257,263],      L"v_r  \;[\mathrm{km s}^{-1}]"],
   #[:ω_f,     degree,   [0,360],        L"\omega_f  \;[\mathrm{rad}]"],
])

names, scaling, ranges, labels = plotting_props_obs_check

##

set_theme!(CornerPlotting.default_theme())

corner_plot = CornerPlotting.CornerPlot(results,[:P_f, :e_f])#, Symbol("P"), Symbol("e")])








##
corner_plot.fig

##

##

@profile create_corner_plot(results, plotting_props_obs_check,
#f = create_corner_plot(results, plotting_props_obs_check,
    supertitle=system_name
    *" - observables",
    )

##

f = create_corner_plot(results, plotting_props_obs_check,
    supertitle=system_name
    *" - observables",
    )


##
save(system_fname*"_observables.png", f)

f

##
plotting_props_obs_check = SideKicks.createPlottingProps([
   #[:m1_f,    m_sun,    [15,40],        L"M_1\;[M_{\odot}]"],
   #[:m2_f,    m_sun,    [0,25],         L"M_2  \;[M_{\odot}]"],
   #[:P_f,     day,      [10.3,10.45],   L"P_f\;[\mathrm{days}]"],
   #[:e_f,     1,        [0,0.1],        L"e_f"],
   #[:K1,      km_per_s, [77,90],        L"K_1  \;[\mathrm{km s}^{-1}]"],
   [:v_N,     km_per_s,  [20, 160],        L"v_N  \;[\mathrm{km s}^{-1}]"],
   [:v_E,     km_per_s,  [20, 160],        L"v_E  \;[\mathrm{km s}^{-1}]"],
   [:v_r,     km_per_s,  [20, 160],        L"v_r  \;[\mathrm{km s}^{-1}]"],
   [:ω_f,     degree,    [0, 360],        L"\omega_f  \;[\mathrm{deg}]"],
])

f = create_corner_plot(results, plotting_props_obs_check,
    supertitle=system_name
    *" - observables 2",
    )
save(system_fname*"_observables2.png", f)

f

##

plotting_props = SideKicks.createPlottingProps([
    [:m2_f,  m_sun,    [0, 20],   L"M_2  \;[M_{\odot}]"],
    [:dm2,   m_sun,    [0, 15],  L"ΔM_2  \;[M_{\odot}]"],
    [:P_f,   day,      [49.9, 50.1],    L"P_f\;[\mathrm{days}]"],
    [:vkick, km_per_s, [0, 150],   L"v_{kick}  \;[\mathrm{km s}^{-1}]"],
    [:vsys,  km_per_s, [0, 100],  L"v_{\mathrm{sys}} \;[\mathrm{km s}^{-1}]"],
])

f = create_corner_plot(results, plotting_props,
    supertitle=system_name
    *" - derived quantities",
    )

save(system_fname*"_derived.png", f)

f
