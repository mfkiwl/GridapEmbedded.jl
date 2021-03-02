using DrWatson
@quickactivate "GridapEmbedded"

using Makie
using AbstractPlotting
using AbstractPlotting.MakieLayout

using DataFrames

function get_raw_data(k::Int,d::Int,t::Int,s::Int,g::Int)
  case = ( (df.k .== k) .& (df.d .== d) .& (df.t .== t) .& (df.s .== s) .& (df.g .== g) )
  ndofs = df[case,:udofs]
  el2Ms = df[case,:el2M]
  eh1Ms = df[case,:eh1M]
  kopMs = df[case,:kopM]
  iteMs = df[case,:iteM]
  el2Ns = df[case,:el2N]
  eh1Ns = df[case,:eh1N]
  kopNs = df[case,:kopN]
  iteNs = df[case,:iteN]
  ndofs, el2Ms, eh1Ms, kopMs, iteMs, el2Ns, eh1Ns, kopNs, iteNs
end

function plot(case::Dict)

  outer_padding = 30
  scene, layout = layoutscene(outer_padding, resolution = (1000, 400),
                              backgroundcolor = RGBf0(0.99, 0.99, 0.99))

  axL = layout[1,1] = LAxis(scene)
  axH = layout[1,2] = LAxis(scene)
  axC = layout[1,3] = LAxis(scene)
  axI = layout[1,4] = LAxis(scene)

  @unpack d,t,s,g = case

  k = 1
  ndofs, el2Ms, eh1Ms, kopMs, iteMs, el2Ns, eh1Ns, kopNs, iteNs = get_raw_data(k,d,t,s,g)

  data = (log10.(ndofs),log10.(el2Ms))
  lineML21 = lines!(axL,data,color=:red,linewidth=2)
  scatML21 = scatter!(axL,data,color=:red,marker=:circle,markersize=6.0)

  data = (log10.(ndofs),log10.(eh1Ms))
  lineMH11 = lines!(axH,data,color=:red,linewidth=2)
  scatMH11 = scatter!(axH,data,color=:red,marker=:xcross,markersize=6.0)

  data = (log10.(ndofs),log10.(kopMs))
  lineMCN1 = lines!(axC,data,color=:red,linewidth=2)
  scatMCN1 = scatter!(axC,data,color=:red,marker=:diamond,markersize=6.0)

  data = (log10.(ndofs),log10.(iteMs))
  lineMIT1 = lines!(axI,data,color=:red,linewidth=2)
  scatMIT1 = scatter!(axI,data,color=:red,marker=:utriangle,markersize=6.0)

  data = (log10.(ndofs),log10.(el2Ns))
  lineNL21 = lines!(axL,data,color=:blue,linestyle=:dash,linewidth=2)
  scatNL21 = scatter!(axL,data,color=:blue,marker=:circle,markersize=6.0)

  data = (log10.(ndofs),log10.(eh1Ns))
  lineNH11 = lines!(axH,data,color=:blue,linestyle=:dash,linewidth=2)
  scatNH11 = scatter!(axH,data,color=:blue,marker=:xcross,markersize=6.0)

  data = (log10.(ndofs),log10.(kopNs))
  lineNCN1 = lines!(axC,data,color=:blue,linestyle=:dash,linewidth=2)
  scatNCN1 = scatter!(axC,data,color=:blue,marker=:diamond,markersize=6.0)

  data = (log10.(ndofs),log10.(iteNs))
  lineNIT1 = lines!(axI,data,color=:blue,linestyle=:dash,linewidth=2)
  scatNIT1 = scatter!(axI,data,color=:blue,marker=:utriangle,markersize=6.0)

  k = 2
  ndofs, el2Ms, eh1Ms, kopMs, iteMs, el2Ns, eh1Ns, kopNs, iteNs = get_raw_data(k,d,t,s,g)

  data = (log10.(ndofs),log10.(el2Ms))
  lineML22 = lines!(axL,data,color=:red,linewidth=2)
  scatML22 = scatter!(axL,data,color=:red,marker=:circle,markersize=6.0)

  data = (log10.(ndofs),log10.(eh1Ms))
  lineMH12 = lines!(axH,data,color=:red,linewidth=2)
  scatMH12 = scatter!(axH,data,color=:red,marker=:xcross,markersize=6.0)

  data = (log10.(ndofs),log10.(kopMs))
  lineMCN2 = lines!(axC,data,color=:red,linewidth=2)
  scatMCN2 = scatter!(axC,data,color=:red,marker=:diamond,markersize=6.0)

  data = (log10.(ndofs),log10.(iteMs))
  lineMIT2 = lines!(axI,data,color=:red,linewidth=2)
  scatMIT2 = scatter!(axI,data,color=:red,marker=:utriangle,markersize=6.0)

  data = (log10.(ndofs),log10.(el2Ns))
  lineNL22 = lines!(axL,data,color=:blue,linestyle=:dash,linewidth=2)
  scatNL22 = scatter!(axL,data,color=:blue,marker=:circle,markersize=6.0)

  data = (log10.(ndofs),log10.(eh1Ns))
  lineNH12 = lines!(axH,data,color=:blue,linestyle=:dash,linewidth=2)
  scatNH12 = scatter!(axH,data,color=:blue,marker=:xcross,markersize=6.0)

  data = (log10.(ndofs),log10.(kopNs))
  lineNCN2 = lines!(axC,data,color=:blue,linestyle=:dash,linewidth=2)
  scatNCN2 = scatter!(axC,data,color=:blue,marker=:diamond,markersize=6.0)

  data = (log10.(ndofs),log10.(iteNs))
  lineNIT2 = lines!(axI,data,color=:blue,linestyle=:dash,linewidth=2)
  scatNIT2 = scatter!(axI,data,color=:blue,marker=:utriangle,markersize=6.0)

  k = 3
  ndofs, el2Ms, eh1Ms, kopMs, iteMs, el2Ns, eh1Ns, kopNs, iteNs = get_raw_data(k,d,t,s,g)

  data = (log10.(ndofs),log10.(el2Ms))
  lineML23 = lines!(axL,data,color=:red,linewidth=2)
  scatML23 = scatter!(axL,data,color=:red,marker=:circle,markersize=6.0)

  data = (log10.(ndofs),log10.(eh1Ms))
  lineMH13 = lines!(axH,data,color=:red,linewidth=2)
  scatMH13 = scatter!(axH,data,color=:red,marker=:xcross,markersize=6.0)

  data = (log10.(ndofs),log10.(kopMs))
  lineMCN3 = lines!(axC,data,color=:red,linewidth=2)
  scatMCN3 = scatter!(axC,data,color=:red,marker=:diamond,markersize=6.0)

  data = (log10.(ndofs),log10.(iteMs))
  lineMIT3 = lines!(axI,data,color=:red,linewidth=2)
  scatMIT3 = scatter!(axI,data,color=:red,marker=:utriangle,markersize=6.0)

  data = (log10.(ndofs),log10.(el2Ns))
  lineNL23 = lines!(axL,data,color=:blue,linestyle=:dash,linewidth=2)
  scatNL23 = scatter!(axL,data,color=:blue,marker=:circle,markersize=6.0)

  data = (log10.(ndofs),log10.(eh1Ns))
  lineNH13 = lines!(axH,data,color=:blue,linestyle=:dash,linewidth=2)
  scatNH13 = scatter!(axH,data,color=:blue,marker=:xcross,markersize=6.0)

  data = (log10.(ndofs),log10.(kopNs))
  lineNCN3 = lines!(axC,data,color=:blue,linestyle=:dash,linewidth=2)
  scatNCN3 = scatter!(axC,data,color=:blue,marker=:diamond,markersize=6.0)

  data = (log10.(ndofs),log10.(iteNs))
  lineNIT3 = lines!(axI,data,color=:blue,linestyle=:dash,linewidth=2)
  scatNIT3 = scatter!(axI,data,color=:blue,marker=:utriangle,markersize=6.0)

  k = 4
  ndofs, el2Ms, eh1Ms, kopMs, iteMs, el2Ns, eh1Ns, kopNs, iteNs = get_raw_data(k,d,t,s,g)

  data = (log10.(ndofs),log10.(el2Ms))
  lineML24 = lines!(axL,data,color=:red,linewidth=2)
  scatML24 = scatter!(axL,data,color=:red,marker=:circle,markersize=6.0)

  data = (log10.(ndofs),log10.(eh1Ms))
  lineMH14 = lines!(axH,data,color=:red,linewidth=2)
  scatMH14 = scatter!(axH,data,color=:red,marker=:xcross,markersize=6.0)

  data = (log10.(ndofs),log10.(kopMs))
  lineMCN4 = lines!(axC,data,color=:red,linewidth=2)
  scatMCN4 = scatter!(axC,data,color=:red,marker=:diamond,markersize=6.0)

  data = (log10.(ndofs),log10.(iteMs))
  lineMIT4 = lines!(axI,data,color=:red,linewidth=2)
  scatMIT4 = scatter!(axI,data,color=:red,marker=:utriangle,markersize=6.0)

  data = (log10.(ndofs),log10.(el2Ns))
  lineNL24 = lines!(axL,data,color=:blue,linestyle=:dash,linewidth=2)
  scatNL24 = scatter!(axL,data,color=:blue,marker=:circle,markersize=6.0)

  data = (log10.(ndofs),log10.(eh1Ns))
  lineNH14 = lines!(axH,data,color=:blue,linestyle=:dash,linewidth=2)
  scatNH14 = scatter!(axH,data,color=:blue,marker=:xcross,markersize=6.0)

  data = (log10.(ndofs),log10.(kopNs))
  lineNCN4 = lines!(axC,data,color=:blue,linestyle=:dash,linewidth=2)
  scatNCN4 = scatter!(axC,data,color=:blue,marker=:diamond,markersize=6.0)

  data = (log10.(ndofs),log10.(iteNs))
  lineNIT4 = lines!(axI,data,color=:blue,linestyle=:dash,linewidth=2)
  scatNIT4 = scatter!(axI,data,color=:blue,marker=:utriangle,markersize=6.0)

  k = 5
  ndofs, el2Ms, eh1Ms, kopMs, iteMs, el2Ns, eh1Ns, kopNs, iteNs = get_raw_data(k,d,t,s,g)

  data = (log10.(ndofs),log10.(el2Ms))
  lineML25 = lines!(axL,data,color=:red,linewidth=2)
  scatML25 = scatter!(axL,data,color=:red,marker=:circle,markersize=6.0)

  data = (log10.(ndofs),log10.(eh1Ms))
  lineMH15 = lines!(axH,data,color=:red,linewidth=2)
  scatMH15 = scatter!(axH,data,color=:red,marker=:xcross,markersize=6.0)

  data = (log10.(ndofs),log10.(kopMs))
  lineMCN5 = lines!(axC,data,color=:red,linewidth=2)
  scatMCN5 = scatter!(axC,data,color=:red,marker=:diamond,markersize=6.0)

  data = (log10.(ndofs),log10.(iteMs))
  lineMIT5 = lines!(axI,data,color=:red,linewidth=2)
  scatMIT5 = scatter!(axI,data,color=:red,marker=:utriangle,markersize=6.0)

  data = (log10.(ndofs),log10.(el2Ns))
  lineNL25 = lines!(axL,data,color=:blue,linestyle=:dash,linewidth=2)
  scatNL25 = scatter!(axL,data,color=:blue,marker=:circle,markersize=6.0)

  data = (log10.(ndofs),log10.(eh1Ns))
  lineNH15 = lines!(axH,data,color=:blue,linestyle=:dash,linewidth=2)
  scatNH15 = scatter!(axH,data,color=:blue,marker=:xcross,markersize=6.0)

  data = (log10.(ndofs),log10.(kopNs))
  lineNCN5 = lines!(axC,data,color=:blue,linestyle=:dash,linewidth=2)
  scatNCN5 = scatter!(axC,data,color=:blue,marker=:diamond,markersize=6.0)

  data = (log10.(ndofs),log10.(iteNs))
  lineNIT5 = lines!(axI,data,color=:blue,linestyle=:dash,linewidth=2)
  scatNIT5 = scatter!(axI,data,color=:blue,marker=:utriangle,markersize=6.0)

  k = 6
  ndofs, el2Ms, eh1Ms, kopMs, iteMs, el2Ns, eh1Ns, kopNs, iteNs = get_raw_data(k,d,t,s,g)

  data = (log10.(ndofs),log10.(el2Ms))
  lineML26 = lines!(axL,data,color=:red,linewidth=2)
  scatML26 = scatter!(axL,data,color=:red,marker=:circle,markersize=6.0)

  data = (log10.(ndofs),log10.(eh1Ms))
  lineMH16 = lines!(axH,data,color=:red,linewidth=2)
  scatMH16 = scatter!(axH,data,color=:red,marker=:xcross,markersize=6.0)

  data = (log10.(ndofs),log10.(kopMs))
  lineMCN6 = lines!(axC,data,color=:red,linewidth=2)
  scatMCN6 = scatter!(axC,data,color=:red,marker=:diamond,markersize=6.0)

  data = (log10.(ndofs),log10.(iteMs))
  lineMIT6 = lines!(axI,data,color=:red,linewidth=2)
  scatMIT6 = scatter!(axI,data,color=:red,marker=:utriangle,markersize=6.0)

  data = (log10.(ndofs),log10.(el2Ns))
  lineNL26 = lines!(axL,data,color=:blue,linestyle=:dash,linewidth=2)
  scatNL26 = scatter!(axL,data,color=:blue,marker=:circle,markersize=6.0)

  data = (log10.(ndofs),log10.(eh1Ns))
  lineNH16 = lines!(axH,data,color=:blue,linestyle=:dash,linewidth=2)
  scatNH16 = scatter!(axH,data,color=:blue,marker=:xcross,markersize=6.0)

  data = (log10.(ndofs),log10.(kopNs))
  lineNCN6 = lines!(axC,data,color=:blue,linestyle=:dash,linewidth=2)
  scatNCN6 = scatter!(axC,data,color=:blue,marker=:diamond,markersize=6.0)

  data = (log10.(ndofs),log10.(iteNs))
  lineNIT6 = lines!(axI,data,color=:blue,linestyle=:dash,linewidth=2)
  scatNIT6 = scatter!(axI,data,color=:blue,marker=:utriangle,markersize=6.0)

  limits!(axL,0.4,3.1,-18.5,0.5)
  axL.xticks = 0.5:0.5:3.0
  axL.yticks = -18.0:2.0:0.0
  axL.xlabel="log10(N^(1/D))"
  axL.ylabel="log10(L2 Abs error)"

  limits!(axH,0.4,3.1,-18.5,0.5)
  axH.xticks = 0.5:0.5:3.0
  axH.yticks = -18.0:2.0:0.0
  axH.xlabel="log10(N^(1/D))"
  axH.ylabel="log10(H1 Abs error)"

  limits!(axC,0.4,3.1,-0.5,13.5)
  axC.xticks = 0.5:0.5:3.0
  axC.yticks = 0.0:1.0:13.0
  axC.xlabel="log10(N^(1/D))"
  axC.ylabel="log10(Condition number)"

  limits!(axI,0.4,3.1,-0.25,4.25)
  axI.xticks = 0.5:0.5:3.0
  axI.yticks = 0.0:0.5:4.0
  axI.xlabel="log10(N^(1/D))"
  axI.ylabel="log10(ML Iters)"

  axL.xticksize = 2.0; axL.yticksize = 2.0
  axL.xticklabelsize = 11.0; axL.yticklabelsize = 11.0
  axL.xlabelsize = 11.0; axL.ylabelsize = 11.0

  axH.xticksize = 2.0; axH.yticksize = 2.0
  axH.xticklabelsize = 11.0; axH.yticklabelsize = 11.0
  axH.xlabelsize = 11.0; axH.ylabelsize = 11.0

  axC.xticksize = 2.0; axC.yticksize = 2.0
  axC.xticklabelsize = 11.0; axC.yticklabelsize = 11.0
  axC.xlabelsize = 11.0; axC.ylabelsize = 11.0

  axI.xticksize = 2.0; axI.yticksize = 2.0
  axI.xticklabelsize = 11.0; axI.yticklabelsize = 11.0
  axI.xlabelsize = 11.0; axI.ylabelsize = 11.0

  mark1 = MarkerElement(color=:gray,marker=:circle,strokecolor=:black)
  mark2 = MarkerElement(color=:gray,marker=:xcross,strokecolor=:black)
  mark3 = MarkerElement(color=:gray,marker=:diamond,strokecolor=:black)
  legmarkers = [ lineML21, lineNL21, mark1, mark2, mark3 ]
  legnames = [ "Modal", "Nodal", "L2(e)", "H1(e)", "κ(A)" ]
  leg = LLegend( scene, legmarkers, legnames, orientation = :vertical )
  leg.labelsize = 11.0
  layout[1,5] = leg

  file = savename("modalC0",case,"png")
  safesave(plotsdir(file),scene)

end

df = collect_results(datadir())
sort!(df,:udofs)

function run()
  params = Dict(
    :d => [2],
    :t => [0,1],
    :s => [0,1],
    :g => [0,1]
  )
  dicts = dict_list(params)
  map(plot,dicts)
end

run()
