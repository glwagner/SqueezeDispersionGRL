# analysis functions

export getkappa, getfinalstate, gettimestep, getcumulants, getgrid, getgridspacing, getstrain,
       reconstructproblem, calccumulants

function getkappa(filename)
  file = jldopen(filename, "r+")
  kap = file["params/kap"]
  close(file)
  kap
end

function reconstructproblem(filename)
  jldopen(filename, "r") do file
     nx = file["grid/nx"]
     ny = file["grid/ny"]
     Lx = file["grid/Lx"]
     Ly = file["grid/Ly"]
      x = file["grid/X"][:, 1]
      y = file["grid/Y"][1, :]
    kap = file["params/kap"]
    eta = file["params/eta"]
      u = file["params/u"]
      v = file["params/v"]
     dt = file["timestepper/dt"]
  end

  c, t, step = getfinalstate(filename)
  x0 = x[1]
  y0 = y[1]

  # Initialize the problem
  g = FourierFlows.TwoDGrid(nx, Lx, ny, Ly; x0=x0, y0=y0, effort=FFTW.MEASURE)
  prob = TracerAdvDiff.ConstDiffProblem(grid=g, kap=kap, eta=eta, u=u, v=v, dt=dt, steadyflow=true)
  TracerAdvDiff.set_c!(prob, c)
  prob.t = t
  prob.step = step

  prob
end

function getfinalstate(filename)
  jldopen(filename, "r") do file
    steps = parse.(keys(file["timeseries/c"]))
    laststep = steps[end]
    c = file["timeseries/c/$laststep"]
    t = file["timeseries/t/$laststep"]
  end
  c, t, laststep
end

function getgrid(filename)
  jldopen(filename, "r") do file
    X = file["grid/X"]
    Y = file["grid/Y"]
  end
  X, Y
end

function getgridspacing(filename)
  jldopen(filename, "r") do file
    x = file["grid/X"][:, 1]
    y = file["grid/Y"][1, :]
    dx = x[2]-x[1]
    dy = y[2]-y[1]
  end
  dx, dy
end

function getstrain(filename)
  jldopen(filename, "r") do file
    t = file["diags/strain/time"]
    S = file["diags/strain/data"]
  end

  # Velocity gradient components
  αdot = map(x -> x[1, 1], S) # extract \dot α from vector S
  γdot = map(x -> x[1, 2], S) # extract \dot γ from vector S
  βdot = map(x -> x[2, 1], S) # extract \dot β from vector S

  t, αdot, βdot, γdot
end

function calccumulants(prob)
  g = prob.grid

  "Integrate `x` over the grid domain."
  integrate(x; dx=g.dx, dy=g.dy) = dx*dy*sum(x)

  @. prob.vars.ch = prob.state.sol
  A_mul_B!(prob.vars.c, g.irfftplan, prob.vars.ch)

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

  C, x⁺, y⁺, m₁₁, m₁₂, m₂₂
end


function getcumulants(filename)
  jldopen(filename, "r") do
    t = file["diags/strain/time"]
    ξ = file["diags/cumulants/data"]
  end
  
  # Cumulants
    C = map(x -> x[1], ξ)
   x⁺ = map(x -> x[2], ξ)
   y⁺ = map(x -> x[3], ξ)
  m₁₁ = map(x -> x[4], ξ)
  m₁₂ = map(x -> x[5], ξ)
  m₂₂ = map(x -> x[6], ξ)

  t, C, x⁺, y⁺, m₁₁, m₁₂, m₂₂
end
