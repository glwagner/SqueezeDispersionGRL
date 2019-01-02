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
    sigexp = extrapolate(sigitp, NaN)
    profsig = sigexp(depth)

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


