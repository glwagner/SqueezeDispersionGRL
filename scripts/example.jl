using Pkg
Pkg.activate("..")

using SqueezeDispersion

# Generate figure 2 data and plot the squeeze dispersion enhancement
#   - runs numerical simulations for delta=0.1:0.1:0.8
#   - saves both simulation results and also diffusivity enhancement
#   - see src/manuscripttools.jl for function definition
deltas, enhancements, setup = genfig2data() 

# Plot the enhancement
plotfig2enhancement(deltas, enhancements, setup)

# Run convergence tests for selected deltas
#= 
for delta in [0.75, 0.8]
  for ny in [256, 512]

    setupkwargs = Dict(
        :delta => delta,
           :nx => 512,
           :ny => ny,
            :H => 1.0,
            :L => 20.0,
          :kap => 1e-5,
          :dcx => 0.2,
          :dcy => 0.05,
       :nsaves => 0,
    :nmessages => 4,
         :name => "test"
    )

    prob, output, diags, setup = setupproblem(; setupkwargs...)
    runproblem(prob, output, diags, setup; nplots=0, showplot=false, saveplot=false)

    h(x) = 1 - setup.delta*sin(2π*x/setup.L)
    println(squeezediffusivity(h; kap=setup.kap, L=setup.L, n=1000) / setup.kap)
    println(calcenhancement(diags[1], setup))

  end
end
=#
