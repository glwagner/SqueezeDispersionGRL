module Samoan

export 
  data2array, 
  cleantempdata, 
  cleansigmadata, 
  simulsort!, 
  movingavg,

  sigma4, 
  latlondist, 
  bindata, 
  binkappa, 
  closestctd,
  getsact,  
  getkappa, 
  getsgth,

  nanmean, 
  nanmaximum, 
  nanminimum,
  @rmnans!, 
  @nan2missing,
  nan2missing, 
  nan2na!, 
  nanoverturns, 
  nanextremes!, 
  
  unpackctddata, 
  unpacksection9,
  smoothoverturns!, 
  vmplayeranalysis, 
  effectivekappa,
  samoanavg,

  pressenter, 
  @initarrays, 
  @makeflat, 

  ABSMAXSIGMA, 
  ABSMINSIGMA, 
  ABSMAXTEMP, 
  ABSMINTEMP

using 
  Statistics,
  PyPlot, 
  Interpolations, 
  Missings, 
  Printf

const ABSMAXSIGMA = 100.0
const ABSMINSIGMA = 10.0
const ABSMAXTEMP = 100.0
const ABSMINTEMP = -1.0
const ABSMINKAPPA = 0.0
const ABSMAXKAPPA = 100.0

include("samoan_utils.jl")
include("data_wrangling.jl")

samoanavg(Φ) = dropdims(nanmean(Φ, 2), dims=2)

"""
    bindata(xedges, x, d)

Put the data in d(x) into x-bins.

Bins and edges geometery:

     bin     bin     bin     bin
  # ----- # ----- # ----- # ----- #
 edge    edge    edge    edge    edge

Thus nbins = nedges-1.

Returns 
=======

binsizes : the number of data in each bin
dbinned : the 'binned data', which is average of the data in each bin
"""
function bindata(xedges, x, d)

  if size(x) != size(d)
    throw("The size of the coordinates and data must be the same.")
  elseif maximum(d) < minimum(xedges) || minimum(d) > maximum(xedges)
    #@info "There is no data in the range spanned by bins."
  end

  nbins = length(xedges)-1
  dbinned = zeros(eltype(d), nbins)
  binsizes = zeros(Int, nbins)

  for i = 1:nbins
    indices = xedges[i] .<= x .<= xedges[i+1] # get indices of data within bins
    binsizes[i] = sum(indices)
    dbinned[i] = mean(d[indices]) # binned data is the mean
  end

  binsizes, dbinned
end

function smoothoverturns!(sigma)
  for i = 2:length(sigma)
    sigma[i] = (sigma[i] < sigma[i-1]) ? sigma[i-1] : sigma[i]
  end
end

function on_density_grid(u, ρ, ρ_grid)
  goodidx = .!isnan.(ρ)
  u_good = u[goodidx]
  ρ_good = ρ[goodidx]
  bindata(ρ_grid, ρ_good, u_good)
end

function thickness(ρ, depth, ρ_grid)
  goodidx = .!isnan.(ρ)
  ρ_good = ρ[goodidx]
  d_good = depth[goodidx]

  # Smooth overturns and calculate thickness across density layers
  smoothoverturns!(ρ_good)
  ρ_interp = interpolate((ρ_good,), d_good, Gridded(Linear()))
  ρ_extrap = extrapolate(ρ_interp, NaN)
  d_grid = ρ_extrap(ρ_grid)
  d_grid[2:end] - d_grid[1:end-1]
end

function vmp_on_density_grid(vmp, ρ_grid)
  nvmp = length(vmp)
  nρ = length(ρ_grid)-1
  
  εᵇ = zeros(nρ, nvmp)
  hᵇ = zeros(nρ, nvmp)
  samples = zeros(nρ, nvmp)

  i = 1
  for (name, prof) in vmp
    ρ = prof["sgth4"]
    d = prof["depth"]
    ε = prof["depth"]

    samples[:, i], εᵇ[:, i] = on_density_grid(ε, ρ, ρ_grid)
    hᵇ[:, i] = thickness(ρ, d, ρ_grid)

    i += 1
  end

  samples, εᵇ, hᵇ
end

"""
    effectivekappa(Γ, ε, h)

Returns κₑ = <h><Γε>, where <Φ> denotes a horizontal average
over the second dimension of Φ (an averge along the 
length of the Samoan passage).

Args
====
Γ : the mixing coefficient
ε : a 2D array of dissipation in (ρ, x) space, where x is distance along
    the Samoan passage and ρ is density.
h : a 2D array of thicknesses.
"""
function effectivekappa(Γ, ε, h; avgh=samoanavg(h))
  flux = @. Γ*ε
  avgh .* samoanavg(flux)
end

"""
    effectivekappa(κ, h)

Returns κₑ = <h><κ/h>, where <Φ> denotes a horizontal average
over the second dimension of Φ (an averge along the 
length of the samoan passage).

Args
====
κ : a 2D array of diffusivity in (ρ, x) space, where x is distance along
    the Samoan passage and ρ is density.
h : a 2D array of thicknesses.
"""
function effectivekappa(κ, h; avgh=samoanavg(h)) 
  flux = @. κ/h
  flux[.!isfinite.(flux)] .= NaN
  avgh .* samoanavg(flux)
end

function localkappa(σ, ε, z; ρ0=1033, dsmooth=20, g=9.81, Γ=0.2)
  σ_avg = movingavg(σ, dsmooth)
  ε_avg = movingavg(ε, dsmooth)
  z_avg = movingavg(z, dsmooth)

  σz = dz(σ_avg, z)
  N² = -g/ρ0 * σz

  return Γ .* ε_avg ./ N² 
end

function movingavg(u, d)
  uavg = zeros(floor(Int, length(u)/d))
  i = j = 1
  i2(k) = k + d - 1
  while i2(i) <= length(u)
    uavg[j] = mean(u[i:i2(i)])
    i += d
    j += 1
  end
  return uavg
end

function vmplayeranalysis(vmp, sigmaedges)
  nvmp = length(vmp)
  nlayers = length(sigmaedges)-1

  layerthickness = zeros(nlayers, nvmp)
  layerkappa = zeros(nlayers, nvmp)
  layerepsilon = zeros(nlayers, nvmp)
  layersamples = zeros(nlayers, nvmp)

  i = 1
  for (name, prof) in vmp
    sigma = prof["sgth4"]
    depth = prof["depth"]
    kappa = prof["kp"]
    epsilon = prof["epsilon"]
    
    # Restrict data to regions with valid density and dissipation
    goodidx = epsilon .> 0
    goodidx .*= .!isnan.(sigma)
    
    goodkappa = kappa[goodidx]
    gooddepth = depth[goodidx]
    goodsigma = sigma[goodidx]
    goodepsilon = epsilon[goodidx]
    
    # Smooth overturns and calculate thickness across density layers
    smoothoverturns!(goodsigma)
    interpsigma = interpolate((goodsigma,), gooddepth, Gridded(Linear()))
    extrapsigma = extrapolate(interpsigma, NaN)
    edgedepths = extrapsigma(sigmaedges)
    layerthickness[:, i] = edgedepths[2:end] - edgedepths[1:end-1]
    
    _, layerkappa[:, i] = bindata(sigmaedges, goodsigma, goodkappa)
    layersamples[:, 1], layerepsilon[:, i] = bindata(sigmaedges, goodsigma, goodepsilon)

    i += 1
  end 

  layersamples, layerkappa, layerepsilon, layerthickness
end

end # module
