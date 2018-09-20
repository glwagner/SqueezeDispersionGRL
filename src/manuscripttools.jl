defaultsetupkwargs = Dict(
   :delta => 0.1,
      :nx => 512,
      :ny => 512,
       :H => 1.0,
       :L => 20.0,
     :kap => 1e-4,
     :dcx => 0.2,
     :dcy => 0.05,
  :nsaves => 10,
    :name => "figure2"
)

"Generate the data used to create figure 2."
function genfig2data(; deltas=0.1:0.1:0.8, showplot=false, saveplot=false, nplots=0, saveenhancements=true)
  enhancements = fill(0.0, size(deltas))
  setupkwargs = deepcopy(defaultsetupkwargs)

  for (i, delta) in enumerate(deltas)
    setupkwargs[:delta] = delta # change delta
    plotdir = joinpath("plots", "figure2", @sprintf("delta%02d", 100delta))

    prob, output, diags, setup = setupproblem(; setupkwargs...)
    runproblem(prob, output, diags, setup; nplots=nplots, showplot=showplot, saveplot=saveplot, plotdir=plotdir,
               savedata=true)

    h(x) = 1 - setup.delta*sin(2π*x/setup.L)
    println(squeezediffusivity(h; kap=setup.kap, L=setup.L, n=1000) / setup.kap)
    println(calcenhancement(diags[1], setup))

    enhancements[i] = calcenhancement(diags[1], setup)
  end

  if saveenhancements
    resultfilename = joinpath("data", "figure2enhancementdata.jld2")
    @save resultfilename deltas enhancements
  end

  defaultsetup = SinusoidalBarotropicSetup(; defaultsetupkwargs...)

  deltas, enhancements, defaultsetup
end


"Plot the enhancement associated with squeeze disperison in the simulation and theory."
function plotfig2enhancement(numericaldeltas, numericalenhancements, basesetup)
  L = basesetup.L
  kap = basesetup.kap

  theorydeltas = 0.0:0.01:0.8
  theoryenhancements = fill(0.0, size(theorydeltas))
  for (i, del) in enumerate(theorydeltas)
    h(x) = 1 - del*sin(2π*x/L)
    theoryenhancements[i] = constantkappaenhancement(del)
  end

  hotspotenhancements = @. 1/(1-theorydeltas)

  fig = figure()
  plot(numericaldeltas, numericalenhancements, "o"; color="C0", markerfacecolor="none", linewidth=4, 
       label="Simulation")
  plot(theorydeltas, theoryenhancements, "-"; color="C0", linewidth=2, label="Theory")
  plot(theorydeltas, hotspotenhancements, "-"; color="C1", linewidth=2, label="Hotspot")
  xlabel(L"Relative bathymetric height, $\Delta$")
  ylabel(L"Enhancement $\kappa_s / \kappa$")
  show()

  fig
end
