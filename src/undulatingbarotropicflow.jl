"Parameters for a tracer advection-diffusion problem over sinusoidal bathymetry."
@with_kw struct SinusoidalBarotropicSetup
  nx = 512
  ny = 256
  delta = 0.1
  H = 1.0
  L = 10.0
  U = 1.0
  dcx = 0.2
  dcy = 0.05
  kap = 1e-5
  eta = 1e-5
  totsteps = 10
  intsteps = 10
  nsaves = 0
  nmessages = 0
  name = "test"
end

function setupmsg(setup)
  msg = "\nSinusodal barotropic problem \"$(setup.name)\" with \n"
  msg *= @sprintf(" delta: %0.4f\n", setup.delta)
  msg *= @sprintf("nx, ny: %d, %d\n", setup.nx, setup.ny)
  msg *= @sprintf("   kap: %.2e\n", setup.kap)
  msg *= @sprintf("   eta: %.2e\n", setup.eta)
  msg *= @sprintf("    dcx: %.2e\n", setup.dcx)
  msg *= @sprintf("    dcy: %.2e\n", setup.dcy)
  msg
end

getfilename(name, nx, ny, delta) = @sprintf("%s_nx%04d_ny%04d_delta%02d.jld2", name, nx, ny, 100delta)

"Construct a TracerAdvDiff problem for a tracer patch advected by barotropic flow over sinusoidal bathymetry."
function setupproblem(;
      delta = 0.1,        # Relative bathymetric height
         nx = 512,        # horizontal resolution 
         ny = 256,        # vertical resolution 
          H = 1.0,        # domain height
          L = 10.0,       # domain width
          U = 1.0,        # barotropic transport
        dcx = 0.2,        # initial tracer width
        dcy = 0.05,       # initial tracer height
        kap = 1e-5,       # vertical diffusivity
        eta = kap,        # horizontal diffusivity
        CFL = 0.2,        # CFL condition for determining time-step
      dmask = 0.01*H,     # width of velocity mask at channel walls
  nmessages = 10,
     nsaves = 10,        
timestepper = "RK4",
       name = "test",
    savedir = joinpath("data")
  )

  if !isdir(savedir); mkpath(savedir); end

  # Get filename, path, and remove output if exists
  filename = getfilename(name, nx, ny, delta)
  filepath = joinpath(savedir, filename)
  if isfile(filepath); rm(filepath); end # Remove output if it exists

  # Pad domain horizontally and vertically
  Lx = 3L/2        
  Ly = H+delta+4dmask 

  # Get the total number of steps and the timestep
  dx = Lx/nx           
  dy = Ly/ny          
  umax = U / (H*(1-delta)) 
  vmax = delta*U/H        
  dt = minimum([CFL*dx/umax, CFL*dy/vmax]) # first guess at timestep
  tf = H*L/U                                # final time
  totsteps = tf/dt

  # Refine timestep to fit intermediate steps
  if nsaves != 0
    intsteps = ceil(Int, totsteps/nsaves)    # Number of snapshots and diags
    totsteps = intsteps*nsaves
    dt = tf/totsteps
  elseif nmessages != 0
    intsteps = ceil(Int, totsteps/nmessages)
    totsteps = intsteps*nmessages
    dt = tf/totsteps
  else
    totsteps = ceil(Int, totsteps)
    intsteps = totsteps
    dt = tf/totsteps
  end

  # Barotropic velocity field associated with bathymetric channel flow
  u, v = sinusoidalbarotropicflow(delta, H, L, U, dmask)
  
  # Initial condition
  x₀, y₀ = 0.0, -H/2
  ci(x, y) = exp( -(x-x₀)^2/(2*dcx^2) - (y-y₀)^2/(2*dcy^2) ) / (2π*dcx*dcy)

  # Initialize the problem
  g = FourierFlows.TwoDGrid(nx, Lx, ny, Ly; y0=-H-delta-2dmask, x0=-L/4, effort=FFTW.MEASURE)
  prob = TracerAdvDiff.ConstDiffProblem(grid=g, kap=kap, eta=eta, u=u, v=v, dt=dt, stepper=timestepper, 
                                        steadyflow=true) 
  TracerAdvDiff.set_c!(prob, ci)

  # Output
  getc(prob) = irfft(prob.state.sol, prob.grid.nx)
  output = Output(prob, filepath, (:c, getc))

  # Diagnostics
  integrate(x; dx=prob.grid.dx, dy=prob.grid.dy) = dx*dy*sum(x)

  """ 
  Return a vector of the strain tensor calculated at the center of
  mass of the tracer blob.  Note that since u = U/h, we have

    ux = - U hₓ / h^2        vx = y U ( hₓₓ / h^2 - 2*hₓ^2/h^3 )
    uy = 0                   vy = u hₓ / h^2  

  where hₓ = ∂h/∂x and hₓₓ = ∂²h/∂x².
  """
  h, hₓ, hₓₓ = sinusoidalbathymetry(delta, H, L)
  function straintensor(prob)
    g = prob.grid

    @. prob.vars.ch = prob.state.sol
    ldiv!(prob.vars.c, g.rfftplan, prob.vars.ch)

    # Xbar and Ybar
    @. prob.vars.cx = g.X*prob.vars.c
    @. prob.vars.cy = g.Y*prob.vars.c
    C = integrate(prob.vars.c)
    x = integrate(prob.vars.cx) / C
    y = integrate(prob.vars.cy) / C

    # Each strain component is a scalar
    ux = - U*hₓ(x) / h(x)^2
    vx = y*U * ( hₓₓ(x) / h(x)^2 - 2*hₓ(x)^2/h(x)^3 )
    vy = U*hₓ(x) / h(x)^2
    uy = 0.0

    [ux vx; uy vy]
  end

  "Retuns the domain average, first cumulants, and second cumulants in a vector."
  function getcumulants(prob)
    g = prob.grid

    @. prob.vars.ch = prob.state.sol
    ldiv!(prob.vars.c, g.rfftplan, prob.vars.ch)

    # First cumulants
    @. prob.vars.cx = g.X*prob.vars.c
    @. prob.vars.cy = g.Y*prob.vars.c
    C = integrate(prob.vars.c)
    x⁺ = integrate(prob.vars.cx) / C
    y⁺ = integrate(prob.vars.cy) / C

    # 2nd cumulants
    @. prob.vars.cx = (g.X - x⁺)^2 * prob.vars.c
    @. prob.vars.cy = (g.Y - y⁺)^2 * prob.vars.c
    m₁₁ = integrate(prob.vars.cx) / C
    m₂₂ = integrate(prob.vars.cy) / C

    @. prob.vars.cy = (g.X - x⁺) * (g.Y - y⁺) * prob.vars.c
    m₁₂ = integrate(prob.vars.cy) / C

    [C, x⁺, y⁺, m₁₁, m₁₂, m₂₂] 
  end

  cumulantsdiag = Diagnostic(getcumulants, prob; nsteps=totsteps)
     straindiag = Diagnostic(straintensor, prob; nsteps=totsteps)

  diags = [cumulantsdiag, straindiag]
  setup = SinusoidalBarotropicSetup(nx=nx, ny=ny, delta=delta, H=H, L=L, U=U, dcx=dcx, dcy=dcy, kap=kap, 
                                    eta=eta, totsteps=totsteps, intsteps=intsteps, name=name, nsaves=nsaves)

  prob, output, diags, setup
end

"Run the tracer `prob` for `totsteps` with `diags` and `output`, making a plot using `mask` and saving output
every `intsteps`."
function runproblem(prob, output, diags, setup; nplots=10, saveplot=false, showplot=false, savedata=false,
                    plotdir=joinpath("plots"))

  plotfreq = nplots > 0 ? floor(Int, setup.totsteps/nplots) : Inf

  # Diagnostic names
  cumulantsdiag = diags[1]
     straindiag = diags[2]
  
  println(setupmsg(setup))
  println("Running...")

  startwalltime = time_ns()
  saveoutput(output) # save initial condition

  function showmsg() 
    @printf("prog: %.2f, wall: %.2f min, c2: %.2f, tot c: %.2f\n",
            prob.step/setup.totsteps, (time_ns()-startwalltime)*1e-9/60, 
            cumulantsdiag.value[6]/cumulantsdiag.data[1][6],
            cumulantsdiag.value[1]/cumulantsdiag.data[1][1])
    nothing
  end

  showmsg()

  if nplots > 0
    close("all")
    fig, ax = subplots(figsize=(8, 4))
    drawtracerplot(ax, output[:c], prob, setup; name=basename(output.filename)[1:end-5], saveplot=saveplot, 
                   showplot=showplot, plotdir=plotdir)
  end

  while prob.step < setup.totsteps
    stepforward!(prob, diags, setup.intsteps)
    showmsg()

    if savedata
      saveoutput(output)
    end

    if prob.step % plotfreq == 0
      drawtracerplot(ax, output[:c], prob, setup; plotdir=plotdir,
                     name=basename(output.filename)[1:end-5], saveplot=saveplot, showplot=showplot)
    end
  end

  r = calcenhancement(cumulantsdiag, setup)
  println("enhancement: $r")

  savediagnostic(cumulantsdiag, "cumulants", output.filename)
  savediagnostic(straindiag, "strain", output.filename)

  nothing
end

function runproblem(; setupkwargs=Dict{Symbol,Any}(), runkwargs=Dict{Symbol,Any}())
  prob, output, diags, setup = setupproblem(setupkwargs...)
  runproblem(prob, output, diags, setup; runkwargs...)
end

function rundefaultproblem(delta=0.1; showplot=true, kwargs...)
  prob, output, diags, setup = setupproblem(delta=delta)
  runproblem(prob, output, diags, setup; showplot=showplot, kwargs...)
end

"Make a rudimentary plot of the tracer distribution `prob.c`, masked by `mask' and save the result in `plotdir`."
function drawtracerplot(ax, c, prob, setup; name="c", saveplot=false, showplot=false, plotdir=joinpath("plots")) 

  if !isdir(plotdir); mkpath(plotdir); end

  X, Y, step = prob.grid.X, prob.grid.Y, prob.step
  U, H, L, delta = setup.U, setup.H, setup.L, setup.delta
  margin = 0.05

  # Draw tracer
  sca(ax)
  cla()
  pcolormesh(X/L, Y/H, c, cmap="Purples", vmin=0.0)

  # Draw bottom topography
  mask = sinusoidalbathymetrymask(delta, H, L, prob.grid)
  bottom = fill(1.0, size(X))
  bottom[.!mask] .= NaN
  pcolormesh(X/L, Y/H, bottom, vmin=0, vmax=4, cmap="Greys")

  # Draw streamlines
  k₁ = 2π/L      
  hh = @. H*(1 - delta*sin(k₁*X))
  ψ = Y*U./hh
  dψ = 0.1
  contour(X/L, Y/H, ψ, levels=(-1+dψ):dψ:-dψ, linewidths=0.5, colors="k", alpha=0.1, linestyles="solid")

  xlabel(L"x/L")
  ylabel(L"y/H")
  xlim(-margin, 1+margin)
  ylim(-(1+delta*1.1), 0.0)
  aspectratio(H/L)
  
  if saveplot
    savename = @sprintf("%s_%06d.png", name, step)
    savefig(joinpath(plotdir, savename), dpi=240)
  end

  if showplot; pause(1.0); end

  nothing
end

"""
Returns `u`, `v` functions that define barotropic flow over a sinusoidal bottom with
relative bathymetric height `delta`, average depth `H`, and length `L`.
"""
function sinusoidalbarotropicflow(delta, H, L, U, dmask) 
  k₁ = 2π/L      

   h(x) = H*(1 - delta*sin(k₁*x)) # Topography
  hₓ(x) = -H*k₁*delta*cos(k₁*x)   # Topographic x-gradient

  # mask(x, y) goes to zero beyond the solid channel walls
  uppermask(x, y) = tanhmask(y-dmask, dmask) 
  lowermask(x, y) = tanhmask(y+h(x)+dmask, dmask; maskbelow=true)
  mask(x, y) = uppermask(x, y)*lowermask(x, y)

  u(x, y) = U/h(x) #* mask(x, y)
  v(x, y) = y*U*hₓ(x)/h(x)^2 * mask(x, y)

  u, v
end

function sinusoidalbathymetry(delta, H, L)
  k₁ = 2π/L

  # Sinusoidal topo and its derivatives
     h(x) = H*(1 - delta*sin(k₁*x))      
    hₓ(x) = -H*k₁*delta*cos(k₁*x)       
   hₓₓ(x) = H*k₁^2*delta*sin(k₁*x)     

  h, hₓ, hₓₓ 
end

"Returns an array of booleans that masks sinusoidal topography."
function sinusoidalbathymetrymask(delta, H, L, X, Y)
  k₁ = 2π/L      
  h(x) = H*(1 - delta*sin(k₁*x))       # Topography

  nx, ny = size(X)

  mask = Array{Bool}(undef, nx, ny)
  for j=1:ny, i=1:nx;
    mask[i, j] = (Y[i, j] < -h(X[i, j]) || Y[i, j] > 0) ? true : false
  end

  mask
end

sinusoidalbathymetrymask(delta, H, L, g::FourierFlows.TwoDGrid) = sinusoidalbathymetrymask(delta, H, L, g.X, g.Y)

function constantkappaenhancement(delta)
  h(x) = 1 - delta*sin(2π*x)
  squeezediffusivity(h; kap=1, L=1, n=1000)
end
