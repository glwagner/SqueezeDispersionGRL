"Integrate x on a grid with spacing dx, assuming that x is periodic."
∫(x, dx=1) = dx*sum(x)

"""Calculate the squeeze diffusivity associated with `h(x)` over the interval `x=0` to `x=L` using `n` points, 
assuming `h(x)` is periodic."""
function squeezediffusivity(h; kap=1, L=1, n=1000)
  dx = L/n
  x = dx:dx:L 
  # squeeze diffusivity = <h><kap/h>
  ∫(h.(x), dx) * ∫(kap ./ h.(x), dx) / L^2
end

function getkappa(filename)
  file = jldopen(filename, "r+")
  kap = file["params/kap"]
  close(file)
  kap
end

function getxy(filename)
  X, Y = nothing, nothing
  jldopen(filename, "r") do file
    X = file["grid/X"]
    Y = file["grid/Y"]
  end
  X, Y
end

function getgridspacing(filename)
  dx, dy = nothing, nothing
  jldopen(filename, "r") do file
    x = file["grid/X"][:, 1]
    y = file["grid/Y"][1, :]
    dx = x[2]-x[1]
    dy = y[2]-y[1]
  end
  dx, dy
end

function getgridsize(filename)
  X, Y = getxy(filename)
  size(X)
end

function getgridlength(filename)
  X, Y = getxy(filename)
  Lx = X[end, 1] - X[1, 1]
  Ly = Y[end, 1] - Y[1, 1]
  Lx, Ly
end

function getgridorigin(filename)
  X, Y = getxy(filename)
  X[1, 1], Y[1, 1]
end

function getgrid(filename)
  nx, ny = getgridsize(filename) 
  Lx, Ly = getgridlength(filename)
  x0, y0 = getgridorigin(filename)
  FourierFlows.TwoDGrid(nx, Lx, ny, Ly; y0=y0, x0=x0)
end

"Calculate the diffusivity enhancement due to squeeze dispersion."
function calcenhancement(cumulants, setup)
  kapeff = (cumulants.data[end][6] - cumulants.data[1][6]) / (2*(cumulants.time[end]-cumulants.time[1]))
  kapeff / setup.kap
end

"Returns the value of a tanh masking function centered at x₀ and with width Δx."
function tanhmask(x, Δx; maskbelow=false) 
  if maskbelow
    return 0.5*(tanh(x/Δx) + 1)
  else
    return 0.5*(tanh(-x/Δx) + 1)
  end
end
