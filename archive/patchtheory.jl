# Define a barotropic squeezing flow
H = 1.0
L = 2π
U = 1.0

  h(x, ε) = H*(1.0 + ε*sin(2π*x/L))
 hₓ(x, ε) = H*ε*cos(2π*x/L) 
hₓₓ(x, ε) = -H*ε*sin(2π*x/L) 

 u(x, y, ε) = U / h(x, ε)
 w(x, y, ε) = y*U*hₓ(x, ε) / h(x, ε)^2
uₓ(x, y, ε) = -U*hₓ(x, ε) / h(x, ε)^2
wₓ(x, y, ε) = y*U/h(x, ε)^2 * ( hₓₓ(x, ε) - 2*hₓ(x, ε)^2/h(x, ε) )


"""
Use Forward-euler integration to solve for the position ξ, ζ 
and moments m₁₁, m₁₂, m₂₂ of tracer patch in a 2D barotropic flow in x, z.
"""
function patchtheory_barotropicflow(ε, ξ₀, ζ₀, m₁₁₀, m₁₂₀, m₂₂₀, t)

  dt = t[2]-t[1]
  nt = length(t)

  # Initialize
  ξ = zeros(nt)
  ζ = zeros(nt)
  m₁₁ = zeros(nt)
  m₁₂ = zeros(nt)
  m₂₂ = zeros(nt)

  ξ[1] = ξ₀
  ζ[1] = ζ₀
  m₁₁[1] = m₁₁₀
  m₁₂[1] = m₁₂₀
  m₂₂[1] = m₂₂₀

  # βdot = 0.0
  αdot = zeros(nt)
  γdot = zeros(nt)
  
  # Forward Euler integration of particle position and moments.
  for i = 1:nt-1
    # Get gradient components
    αdot[i] = uₓ(ξ[i], ζ[i], ε)
    γdot[i] = wₓ(ξ[i], ζ[i], ε)
    
    # Moment equations. Note that βdot = u_z = 0.
    m₁₁[i+1] = m₁₁[i] + dt*(2κ + 2αdot[i]*m₁₁[i])
    m₁₂[i+1] = m₁₂[i] + dt*γdot[i]*m₁₁[i]
    m₂₂[i+1] = m₂₂[i] + dt*(2κ + 2γdot[i]*m₁₂[i] - 2αdot[i]*m₂₂[i])

    # Update particle position
    ξ[i+1] = ξ[i] + dt*u(ξ[i], ε)
    ζ[i+1] = ζ[i] + dt*w(ξ[i], ζ[i], ε)
  end

  # Get velocity gradient components

  ξ, ζ, m₁₁, m₁₂, m₂₂ 
end
