module Samoan

using MAT, PyPlot, Interpolations, Missings, Printf #, DataArrays

export cleantempdata, cleansigmadata, getsact, getkappa, nan2na!, data2array, getsgth
export simulsort!, pressenter, nan2missing, nanoverturns, @rmnans!, bindata
export @initarrays, @makeflat, @nan2missing, binkappa, closestctd

export latlondist, sigma4, unpackctddata, nanextremes!, smoothoverturns!, unpacksection9,
       vmplayeranalysis, effectivekappa

export nanmean, nanmaximum, nanminimum

export ABSMAXSIGMA, ABSMINSIGMA, ABSMAXTEMP, ABSMINTEMP

const ABSMAXSIGMA = 100.0
const ABSMINSIGMA = 10.0
const ABSMAXTEMP = 100.0
const ABSMINTEMP = -1.0
const ABSMINKAPPA = 0.0
const ABSMAXKAPPA = 100.0

"Put the data in d(x) into x-bins."
function bindata(xedges, x, d)

  if size(x) != size(d)
    throw("The size of the coordinates and data must be the same.")
  elseif maximum(d) < minimum(xedges) || minimum(d) > maximum(xedges)
    warn("There is no data in the range spanned by bins.")
  elseif minimum(d) < minimum(xedges) || maximum(d) > maximum(xedges)
    warn("Data exists outside the range spanned by the bins.")
  end

  # bins and edges geometery:
  # ----- # ----- # ----- # ----- #
  #  bin        edge
  #
  # d = d(x)

  nbins = length(xedges)-1
  dbinned = zeros(eltype(d), nbins)
  binsizes = zeros(Int, nbins)

  for i = 1:nbins
    indexes = xedges[i] .<= x .<= xedges[i+1] # get indexes of data within bins
    binsizes[i] = sum(indexes)
    dbinned[i] = mean(d[indexes]) # binned data is the mean
  end

  binsizes, dbinned
end

function smoothoverturns!(sigma)
  for i = 2:length(sigma)
    sigma[i] = (sigma[i] < sigma[i-1]) ? sigma[i-1] : sigma[i]
  end
end

"Take the maximum of data while ignoring NaNs. "
nanmaximum(x) = maximum(filter(!isnan, x))
nanmaximum(x, y) = mapslices(nanmaximum, x, y)

"Take the minimum of data while ignoring NaNs. "
nanminimum(x) = minimum(filter(!isnan, x))
nanminimum(x, y) = mapslices(nanminimum, x, y)

nanmean(x) = mean(filter(!isnan, x))
nanmean(x, y) = mapslices(nanmean, x, y)

"Calculate the distance between two points at (lat1, lon1) and (lat2, lon2)."
function latlondist(lat1, lon1, lat2, lon2)
  # ref: https://andrew.hedges.name/experiments/haversine/
  R = 6.3178e6 # Earth radius
  dlon = abs(lon1-lon2)
  dlat = abs(lat1-lat2)
  a = sin(dlat/2)^2 + cos(lat1)*cos(lat2)*sin(dlon/2)^2
  c = 2*atan(sqrt(a), sqrt(1-a))
  R*c
end
latlondist(latlon1::Tuple, latlon2::Tuple) = latlondist(latlon1[1], latlon1[2], latlon2[1], latlon2[2])


"""Set the elements of data below lowerlimit and above upperlimit to NaN. Set the corresponding
elements of assocdata to NaN as well."""
function nanextremes!(data...; lowerlimit=-Inf, upperlimit=Inf)
  primarydata = data[1]
  idxhigh = primarydata .> upperlimit
   idxlow = primarydata .< lowerlimit
  for d in data
    d[idxhigh] .= NaN
    d[idxlow] .= NaN
  end
end



function effectivekappa(layerkappa, thickness; pow=1, meanthickness=dropdims(mean(thickness, 2), dims=2))
  kernel = layerkappa ./ thickness.^pow
  kernel[.!isfinite.(kernel)] .= NaN
  dropdims(nanmean(kernel, 2), dims=2) .* meanthickness.^pow
end

function vmplayeranalysis(vmp, sigmaedges)
  nvmp = length(vmp)
  nlayers = length(sigmaedges)-1

  thickness = zeros(nlayers, nvmp)
  layerkappa = zeros(nlayers, nvmp)
  layersamples = zeros(nlayers, nvmp)

  i = 1
  for (name, prof) in vmp
    sigma = prof["sgth4"]
    depth = prof["depth"]
    kappa = prof["kp"]
    
    goodidx = kappa .> 0
    goodidx .*= .!isnan.(sigma)
    
    goodkappa = kappa[goodidx]
    gooddepth = depth[goodidx]
    goodsigma = sigma[goodidx]
    
    smoothoverturns!(goodsigma)
    interpsigma = interpolate((goodsigma,), gooddepth, Gridded(Linear()))
    edgedepths = interpsigma[sigmaedges]
    thickness[:, i] = edgedepths[2:end] - edgedepths[1:end-1]
    
    layersamples[:, i], layerkappa[:, i] = bindata(sigmaedges, goodsigma, goodkappa)

    i += 1
  end 

  layersamples, layerkappa, thickness
end


function unpacksection9(ctddata)
  idx = ctddata["section"] .== 9
  ignore = ((ctddata["cast"] .!= 49) .*  
            (ctddata["cast"] .!= 60) .*
            (ctddata["cast"] .!= 50) .*
            (ctddata["cast"] .!= 51)
           )
  idx .*= ignore

  nsection9 = sum(idx)
  ndepth = length(ctddata["p"][1])

  cast = zeros(1, nsection9)
   lon = zeros(1, nsection9)
   lat = zeros(1, nsection9)
     p = zeros(ndepth, nsection9)
     t = zeros(ndepth, nsection9)
    sp = zeros(ndepth, nsection9)
    sa = zeros(ndepth, nsection9)
  sig4 = zeros(ndepth, nsection9)

  depth = dropdims(ctddata["z"][1], dims=1)

  isec9 = 1
  for (ictd, count) in enumerate(idx)
    if count
      cast[isec9] = ctddata["cast"][ictd] 
       lon[isec9] = ctddata["lon"][ictd] 
       lat[isec9] = ctddata["lat"][ictd] 
         p[:, isec9] = dropdims(ctddata["p"][ictd], dims=1)
         t[:, isec9] = dropdims(ctddata["t1"][ictd], dims=1)
        sp[:, isec9] = dropdims(ctddata["s1"][ictd], dims=1)
      sig4[:, isec9] = dropdims(ctddata["sigma41"][ictd], dims=1)
      isec9 += 1
    end
  end 

  # NaN sigma and temp.
  nanextremes!(sig4, p, t, sp; lowerlimit=ABSMINSIGMA, upperlimit=ABSMAXSIGMA)
  nanextremes!(t, sig4, p, sp; lowerlimit=ABSMINTEMP, upperlimit=ABSMAXTEMP)

  # Get absolute salinity
  #sa = GSW.sa_from_sp.(sp, p, lon, lat)

  #cast, lon, lat, p, t, sp, sa, sig4, depth
  cast, lon, lat, p, t, sp, sig4, depth
end


function getprofilessigma!(profiles, ctdlon, ctdlat, ctdsig)

  nprofiles = length(profiles)
  iprof = 1

  for (name, profile) in profiles
    proflon = profile["lon"]
    proflat = profile["lat"]
    profdepth = dropdims(profile["depth"], dims=2)
       profsp = dropdims(profile["s"], dims=2)
       profkp = dropdims(profile["kp"], dims=2)
        proft = cleantempdata(dropdims(profile["t"], dims=2))

    lon[iprof] = proflon
    lat[iprof] = proflat

    ictd, mindist = closestctd(proflon, proflat, ctdlon, ctdlat)
    icolor = iprof/nprofiles

    closestsig = ctdsig[:, ictd]
      cleansig = closestsig.data[.!closestsig.na]
    cleandepth = ctddepth[.!closestsig.na]

    sigitp = interpolate((cleandepth,), cleansig, Gridded(Linear()))
    profsig = sigitp[depth]

    profile["sigma"] = profsig
  end
end

function getsgth(profile)
  sgth = dropdims(profile["sgth"], dims=2)
  idxhigh = sgth .> ABSMAXSIGMA
   idxlow = sgth .< ABSMINSIGMA
  sgth[idxhigh] = NaN
  sgth[idxlow] = NaN
  sgth
end

function getkappa(profile)
  kappa = dropdims(profile["kp"], dims=2)
  idxhigh = kappa .> ABSMAXKAPPA
   idxlow = kappa .< ABSMINKAPPA
  kappa[idxhigh] = NaN
  kappa[idxlow] = NaN
  kappa
end

function binkappa(profile, σedges; usesgth=false)
   sp = dropdims(profile["s"], dims=2)
    t = dropdims(profile["t"], dims=2)
    p = dropdims(profile["p"], dims=2)
  lon = profile["lon"]
  lat = profile["lat"]
    k = dropdims(profile["kp"], dims=2)

  if usesgth
    σ = getsgth(profile)
  else
    σ = sigma4.(sp, t, p, lon, lat) # get density
  end

  # Remove kappa values for which no corresponding density information exists.
  nans = isnan.(σ)
  σ = σ[.!nans]
  k = k[.!nans]

  bindata(σedges, σ, k)
end


function data2array(adata)
  a = zeros(Float64, size(adata))
  for i in eachindex(adata)
    a[i] = isfinite(adata[i]) ? adata[i] : NaN
  end
  a 
end

function cleansigmadata(sigma)
  idxhigh = sigma .> ABSMAXSIGMA
  idxlow = sigma .< ABSMINSIGMA
  sigma[idxhigh] = NaN
  sigma[idxlow] = NaN
  sigma = data(sigma)
  nan2na!(sigma)
  sigma
end

function cleantempdata(temp)
  idxhigh = temp .> ABSMAXTEMP
  idxlow = temp .< ABSMINTEMP
  temp[idxhigh] = NaN
  temp[idxlow] = NaN
  temp = data(temp)
  nan2na!(temp)
  temp
end


""" Isolate the NaNs in the first array argument and remove corresponding elements in all arrays. """
macro rmnans!(arrays...)
  nans = isnan.(arrays[1])
  expr = Expr(:block)
  append!(expr.args, [:( $(esc(a)) = $(esc(a))[.!nans]; ) for a in arrays])
  expr
end

function simulsort!(arrays...)
  ii = sortperm(arrays[1])
  for a in arrays; a = a[ii]; end
  nothing
end

function pressenter()
  println("Press enter to continue")
  readline(STDIN)
  nothing
end

function nan2missing(a)
  amiss = zeros(Union{Float64,Missing}, size(a))
  amiss .= a
  amiss[isnan.(amiss)] = missing
  amiss
end

macro nan2missing!(arrays...)
  expr = Expr(:block)
  append!(expr.args, [:( $(esc(a)) = nan2missing($(esc(a))); ) for a in arrays])
  expr
end

function nanoverturns(pden)
  npden = length(pden)
  # Assume sorted from shallow to deep (which implies that pden increases)
  i = 2
  while i < npden
    if pden[i] < pden[i-1] # overturn found
      itop = i-1
      ibot = i+1
      while pden[ibot] < pden[ibot-1] && ibot < npden # look for top of overturn
        ibot += 1
      end
      pden[itop:ibot-1] .= NaN # NaN-out the overturn
      i = ibot+1
    else
      i += 1
    end
  end
  if pden[end] < pden[end-1]; pden[end] = NaN; end # deal with bottom cell
  pden
end

function nan2na!(a)
  for i in eachindex(a)
    a[i] = ismissing(a[i]) ? NA : (isna(a[i]) ? NA : (isfinite(a[i]) ? a[i] : NA))
  end
  nothing
end

macro initarrays(vars...)
  expr = Expr(:block)
  append!(expr.args, [:($(esc(var)) = []; ) for var in vars])
  expr
end

macro makeflat(arrays...)
  expr = Expr(:block)
  for a in arrays
    flata = Symbol(:flat, a)
    append!(expr.args, [:($(esc(flata)) = vcat($(esc(a))...); )])
  end
  expr
end

function closestctd(lon, lat, ctdlons, ctdlats)
  # ref: https://andrew.hedges.name/experiments/haversine/
  R = 6.3178e6 # Earth radius
  mind = Inf
  idx = 0
  for i = 1:length(ctdlons)
    ctdlon = ctdlons[i]
    ctdlat = ctdlats[i]
    dlon = abs(lon-ctdlon)
    dlat = abs(lat-ctdlat)
    a = sin(dlat/2)^2 + cos(lat)*cos(ctdlat)*sin(dlon/2)^2
    c = 2*atan2(sqrt(a), sqrt(1-a))
    d = R*c
    if d < mind
      idx = i
      mind = d
    end
  end
  idx, mind
end

#=
import GibbsSeaWater; GSW = GibbsSeaWater

""" Get absolute salinity and conservative temperature. """
function getsact(sp, p, lon, lat)
  sa = zeros(sp)
  ct = zeros(sp)
  for i = 1:length(sp)
    sa[i] = GSW.sa_from_sp(sp[i], p[i], lon, lat)
    ct[i] = GSW.ct_from_t(sa[i], t[i], p[i])
  end
  sa, ct
end

""" Get potential density referenced to 4000 m. """
function sigma4(sp, t, p, lon, lat)
  sa = GSW.sa_from_sp(sp, p, lon, lat)
  ct = GSW.ct_from_t(sa, t, p)
  GSW.sigma4(sa, ct)
end

function unpackctddata(ctddata)
  nprofiles = length(ctddata["lon"])
  ndepth = length(ctddata["p"][1])

  profilelengths = [ length(p) for p in ctddata["p"] ]
  nfullprofiles = length(profilelengths[profilelengths .> 1])

   lon = zeros(1, nfullprofiles)
   lat = zeros(1, nfullprofiles)
     p = zeros(ndepth, nfullprofiles)
     t = zeros(ndepth, nfullprofiles)
    sp = zeros(ndepth, nfullprofiles)
    sa = zeros(ndepth, nfullprofiles)
  sig4 = zeros(ndepth, nfullprofiles)
    
  depth = dropdims(ctddata["z"][1], dims=1)

  ictd = 1
  for i = 1:nprofiles
    if profilelengths[i] > 1
      lon[ictd] = ctddata["lon"][i] 
      lat[ictd] = ctddata["lat"][i] 
         p[:, ictd] = dropdims(ctddata["p"][i], dims=1)
         t[:, ictd] = dropdims(ctddata["t1"][i], dims=1)
        sp[:, ictd] = dropdims(ctddata["s1"][i], dims=1)
      sig4[:, ictd] = dropdims(ctddata["sigma41"][i], dims=1)
      ictd += 1
    end
  end

  # NaN sigma and temp.
  nanextremes!(sig4, p, t, sp; lowerlimit=ABSMINSIGMA, upperlimit=ABSMAXSIGMA)
  nanextremes!(t, sig4, p, sp; lowerlimit=ABSMINTEMP, upperlimit=ABSMAXTEMP)

  # Get absolute salinity
  sa = GSW.sa_from_sp.(sp, p, lon, lat)

  lon, lat, p, t, sp, sa, sig4, depth
end

=#


end # module
