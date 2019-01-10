"Take the maximum of data while ignoring NaNs. "
nanmaximum(x) = maximum(filter(!isnan, x))
nanmaximum(x, y) = mapslices(nanmaximum, x, dims=y)

"Take the minimum of data while ignoring NaNs. "
nanminimum(x) = minimum(filter(!isnan, x))
nanminimum(x, y) = mapslices(nanminimum, x, dims=y)

nanmean(x) = mean(filter(!isnan, x))
nanmean(x, y) = mapslices(nanmean, x, dims=y)

"""
Set the elements of data below lowerlimit and above upperlimit to NaN. 
Set the corresponding elements of assocdata to NaN as well.
"""
function nanextremes!(data...; lowerlimit=-Inf, upperlimit=Inf)
  primarydata = data[1]
  idxhigh = primarydata .> upperlimit
   idxlow = primarydata .< lowerlimit
  for d in data
    d[idxhigh] .= NaN
    d[idxlow] .= NaN
  end
end

"""
Isolate the NaNs in the first array argument and remove corresponding 
elements in all arrays.
"""
macro rmnans!(arrays...)
  expr = Expr(:block)
  append!(expr.args, :(nans = isnan.($(esc(arrays[1])))))
  append!(expr.args, [:( $(esc(a)) = $(esc(a))[.!nans]; ) for a in arrays])
  expr
end

rmnans(a) = a[.!isnan.(a)]

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

function simulsort!(arrays...)
  ii = sortperm(arrays[1])
  for a in arrays; a = a[ii]; end
  nothing
end

function simulsort(arrays...; kwargs...)
  ii = sortperm(arrays[1]; kwargs...)
  (a[ii] for a in arrays)
end

function pressenter()
  println("Press enter to continue")
  readline(STDIN)
  nothing
end

macro makeflat(arrays...)
  expr = Expr(:block)
  for a in arrays
    flata = Symbol(:flat, a)
    append!(expr.args, [:($(esc(flata)) = vcat($(esc(a))...); )])
  end
  expr
end

"Calculate the distance between two points at (lat1, lon1) and (lat2, lon2)."
function latlondist(lat1, lon1, lat2, lon2)
  # ref: https://andrew.hedges.name/experiments/haversine/
  R = 6.3178e6 # Earth radius
  dlon = lon1-lon2
  dlat = lat1-lat2
  a = sin(dlat/2)^2 + cos(lat1)*cos(lat2)*sin(dlon/2)^2
  c = 2*atan(sqrt(a), sqrt(1-a))
  R*c
end

latlondist(latlon1::Tuple, latlon2::Tuple) = latlondist(latlon1[1], latlon1[2], 
                                                        latlon2[1], latlon2[2])

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

"""
    dz!(uz, u, z)

Calculate `uz = du/dz` using centered differences.
"""
function dz!(uz, u, z; rev=true)
  if rev
    uz[1] = (u[1] - u[2]) / (z[1] - z[2])
    uz[end] = (u[end-1] - u[end]) / (z[end-1] - z[end])
    @. @views uz[2:end-1] = (u[1:end-2] - u[3:end]) / (z[1:end-2] - z[3:end])
  else
    uz[1] = (u[2] - u[1]) / (z[2] - z[1])
    uz[end] = (u[end] - u[end-1]) / (z[end] - z[end-1])
    @. @views uz[2:end-1] = (u[3:end] - u[1:end-2]) / (z[3:end] - z[1:end-2])
  end
end

function dz(u, z; kwargs...)
  uz = similar(u)
  dz!(uz, u, z, kwargs...)
  uz
end
