```@meta
CurrentModule = SideKicks
```

# SideKicks.jl

Welcome to the documentation for [SideKicks](https://github.com/orlox/SideKicks.jl), a Julia package for analyzing the observed properties of stellar binaries containing a compact object. The package enables parameter inference on supernova mass loss and natal kicks using advanced statistical methods.

## Overview

SideKicks.jl provides tools for:
- Modeling supernova mass loss and natal kicks in stellar binaries.
- Performing Markov Chain Monte Carlo (MCMC) simulations using Turing.jl.
- Handling observational data and priors for astrophysical systems.
- Visualizing results with corner plots and other diagnostic tools.
- Saving and extracting results in HDF5 format.

## Features

- **Custom Probability Distributions**: Includes `WrappedCauchy` and `ModVonMises` distributions for angular data.
- **Turing Models**: Simplified and general MCMC models for parameter inference.
- **Observational Data Handling**: Tools for managing observed parameters, errors, and units.
- **Visualization**: Generate corner plots for visualizing parameter relationships.
- **Output Storage**: Save and extract results efficiently using HDF5.

## Getting Started

To get started with SideKicks.jl, check out the following sections:
- [Installation](installation.md): Learn how to install the package.
- [Usage](usage.md): Explore examples and workflows for using SideKicks.jl.

## Index

```@index
```

## API Documentation

```@autodocs
Modules = [SideKicks]
```
