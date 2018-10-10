# SqueezeDispersionGRL

This repository contains data, source code, and notebooks for reproducing figures 2 and 3 in the 
manuscript "Squeeze dispersion: modulation of diapycnal mixing by diapycnal strain", submitted to
Geophysical Research Letters. The code is written in Julia and most of the data is stored in the
[JLD2]() data format.

## Data

The observational microstructure and CTD data is stored in `data/samoanpassagedata.jld2`
Bathymetric data for the Samoan passage region is stored in `data/bathy.nc`.
Output from numerical simulations is stored in `data/figure2_kap1e-4_nx0512_ny0512_delta50.jld2` and 
`data/figure2enhancements.jld2`.

## Script to reproduce simulation data

A script to reproduce the simulation data is provided in `scripts/example.jl`.

## Notebooks to produce plots

Jupyter notebooks to produce plots for the paper are provided in `notebooks/`.
The notebook `example.ipynb` produces figure 2. The notebooks `samoanoverview.ipynb` and `samoanresults.ipynb`
produce the top and bottom of figure 3. 

## Author

The code in this repository was written by [Gregory L. Wagner](). The observational data was trimmed from a larger CTD
and microstructure dataset provided by [Glenn S. Carter]() and [Gunnar Voet]().

[JLD2]: https://www.github.com/JLD2.jl
[Gregory L. Wagner]: https://glwagner.github.io
[Glenn S. Carter]:
[Gunnar Voet]:
