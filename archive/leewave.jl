""" 
Returns u, v for a lee wave forming in steady flow with velocity U, 
wavenumbers k, m, and amplitude h, so that

  ψ = U y - U h cos(k x + m y); u = ∂y psi, v = - ∂x psi

  evanescent
    ψ(x, y) = U*y - U*h*cos(k*x)*exp(-m*y)
    u(x, y) = U + m*U*h*cos(k*x)*exp(-m*y)
    v(x, y) = -k*U*h*sin(k*x)*exp(-m*y)

"""
function leewaveflow(U, h, k, m; evanescent=false)

  T = U*h
  
  Φ(x, y) = k*x + m*y
  ψ(x, y) = U*y - T*cos(Φ(x, y))

  u(x, y) = U + m*T*sin(Φ(x, y))
  v(x, y) =   - k*T*sin(Φ(x, y))

  ux(x, y) =  m*k*T*cos(Φ(x, y))
  uy(x, y) =  m^2*T*cos(Φ(x, y))
  vx(x, y) = -k^2*T*cos(Φ(x, y))

  ψ, u, v, ux, uy, vx
end


function leewavepatchrun(U, ε, k, m, δ, ξ₀, ζ₀)

  L = 2π/k
  H = 2π/m
  dt = 0.005*L/U

  ψ, u₁, v₁, ux₁, uy₁, vx₁ = leewaveflow(U, ε, k, m)

  # Set up and solve tracer patch problem
   u★(x, y, t) = u₁(x, y)
   v★(x, y, t) = v₁(x, y)
  ux★(x, y, t) = ux₁(x, y)
  uy★(x, y, t) = uy₁(x, y)
  vx★(x, y, t) = vx₁(x, y)

  prob = TracerPatchEqn.AnalyticFlowProblem(
    κ=κ, η=η, u=u★, v=v★, ux=ux★, uy=uy★, vx=vx★, dt=dt)

  TracerPatchEqn.set_position!(prob, ξ₀, ζ₀)
  TracerPatchEqn.set_moments!(prob, δ^2, 0.0, δ^2)

  ξ₋₁ = prob.vars.sol[1]
  while prob.vars.sol[1] < (ξ₀ + L)
    ξ₋₁ = prob.vars.sol[1]
    stepforward!(prob)
  end

  # Numerical calculation of end time (seems to be T = L/U
  Δ = (prob.vars.sol[1] - (ξ₀ + L))/(prob.vars.sol[1] - ξ₋₁)  # frac overshoot
  T = prob.t - dt*Δ
  dt = T / nsteps

  # Restart
  prob = TracerPatchEqn.AnalyticFlowProblem(
    κ=κ, η=η, u=u★, v=v★, ux=ux★, uy=uy★, vx=vx★, dt=dt)
  TracerPatchEqn.set_position!(prob, ξ₀, ζ₀)
  TracerPatchEqn.set_moments!(prob, δ^2, 0.0, δ^2)

  getsol(prob) = deepcopy(prob.vars.sol)
  diag = Diagnostic(getsol, prob, nsteps=nsteps)
  stepforward!(prob, [diag], nsteps=nsteps)

  # Extract data
  t = diag.time
  ξ = map(x -> x[1], diag.data)
  ζ = map(x -> x[2], diag.data)
  s = map(x -> x[3], diag.data)
  θ = map(x -> x[4], diag.data)
  σ = map(x -> x[5], diag.data)

  t, ξ, ζ, s, θ, σ
end

