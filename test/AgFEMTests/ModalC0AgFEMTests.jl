using Gridap
using Gridap.Fields
using Gridap.Geometry
using GridapEmbedded
using GridapEmbedded.AgFEM
using LinearAlgebra: cond

using Makie
using AbstractPlotting
using AbstractPlotting.MakieLayout

using AlgebraicMultigrid
import IterativeSolvers: cg

u(x) = (x[1] - x[2])^6
f(x) = -Δ(u)(x)
ud(x) = u(x)

const R = 0.42

function run(n::Int,order::Int)

  # geom = disk(R,x0=Fields.Point(0.5,0.5))
  # partition = (n,n)
  # domain = (0,1,0,1)

  geom = sphere(R,x0=Fields.Point(0.5,0.5,0.5))
  partition = (n,n,n)
  domain = (0,1,0,1,0,1)

  bgmodel = CartesianDiscreteModel(domain,partition)
  h = (domain[2]-domain[1])/n
  cutgeo = cut(bgmodel,geom)
  model = DiscreteModel(cutgeo)

  # γ₀ = 0.65
  # eᵧ = (3-2/num_dims(model))/(2*order+1-2/num_dims(model))
  # γ = γ₀^eᵧ
  # strategy = AggregateCutCellsByThreshold(γ)
  strategy = AggregateAllCutCells()

  aggregates = aggregate(strategy,cutgeo,geom)
  bboxes = compute_cell_to_dface_bboxes(model,aggregates)

  Ω_bg = Triangulation(bgmodel)
  Ω = Triangulation(cutgeo)
  Γ = EmbeddedBoundary(cutgeo)

  n_Γ = get_normal_vector(Γ)

  D = num_dims(model)
  cutdeg, degree = 2*order, 2*order
  dΩ = Measure(Ω,cutdeg,degree)
  dΓ = Measure(Γ,cutdeg)

  γd = 2.5*order^2

  a(u,v) =
    ∫( ∇(u)⋅∇(v) ) * dΩ +
    ∫( (γd/h)*u*v  - v*(n_Γ⋅∇(u)) - (n_Γ⋅∇(v))*u ) * dΓ

  l(v) =
    ∫( v*f ) * dΩ +
    ∫( (γd/h)*v*ud - (n_Γ⋅∇(v))*ud ) * dΓ

  l2(u) = sqrt(sum( ∫( u*u )*dΩ ))
  h1(u) = sqrt(sum( ∫( u*u + ∇(u)⋅∇(u) )*dΩ ))

  reffe = ReferenceFE(modalC0,Float64,order,bboxes)
  Vstd = TestFESpace(model,reffe,conformity=:H1)
  V = AgFEMSpace(Vstd,aggregates,modalC0)
  U = TrialFESpace(V)

  op = AffineFEOperator(a,l,U,V)
  kopM = cond(get_matrix(op),1) # https://github.com/JuliaLang/julia/issues/6485

  A = get_matrix(op)
  b = get_vector(op)

  p = aspreconditioner(smoothed_aggregation(A))
  ur, ch = cg(A,b,Pl=p,log=true,maxiter=10000,reltol=1.0e-16)
  iteM = ch.iters

  uh = FEFunction(U,ur)
  e = u - uh
  el2M = l2(e)
  eh1M = h1(e)

  reffe = ReferenceFE(lagrangian,Float64,order)
  Vstd = TestFESpace(model,reffe,conformity=:H1)
  V = AgFEMSpace(Vstd,aggregates,lagrangian)
  U = TrialFESpace(V)

  op = AffineFEOperator(a,l,U,V)
  kopN = cond(get_matrix(op),1) # https://github.com/JuliaLang/julia/issues/6485

  A = get_matrix(op)
  b = get_vector(op)

  p = aspreconditioner(smoothed_aggregation(A))
  ur, ch = cg(A,b,Pl=p,log=true,maxiter=10000,reltol=1.0e-16)
  iteN = ch.iters

  uh = FEFunction(U,ur)
  e = u - uh
  el2N = l2(e)
  eh1N = h1(e)

  num_free_dofs(U)^(1/D), el2M, eh1M, kopM, iteM, el2N, eh1N, kopN, iteN

end

function conv_test(ns,order)

  ndofs = Float64[]
  el2Ms = Float64[]
  eh1Ms = Float64[]
  kopMs = Float64[]
  iteMs = Float64[]
  el2Ns = Float64[]
  eh1Ns = Float64[]
  kopNs = Float64[]
  iteNs = Float64[]

  for n in ns

    @time ndof, el2M, eh1M, kopM, iteM, el2N, eh1N, kopN, iteN = run(n,order)
    println("completed: ",n," ",order)

    push!(ndofs,ndof)
    push!(el2Ms,el2M)
    push!(eh1Ms,eh1M)
    push!(kopMs,kopM)
    push!(iteMs,iteM)
    push!(el2Ns,el2N)
    push!(eh1Ns,eh1N)
    push!(kopNs,kopN)
    push!(iteNs,iteN)

  end

  (ndofs, el2Ms, eh1Ms, kopMs, iteMs, el2Ns, eh1Ns, kopNs, iteNs)

end

function plot()

  outer_padding = 30
  scene, layout = layoutscene(outer_padding, resolution = (1000, 400),
                              backgroundcolor = RGBf0(0.99, 0.99, 0.99))

  axL = layout[1,1] = LAxis(scene)
  axH = layout[1,2] = LAxis(scene)
  axC = layout[1,3] = LAxis(scene)
  axI = layout[1,4] = LAxis(scene)

  # nes = [ [8,16,32,64,128,256,512], [6,12,24,48,96,192], [6,12,24,48,96], [6,12,24,48,96], [6,12,24,48] ]
  nes = [ [8,16,32,64], [6,12,24], [6], [6,12,24], [6,12,24,48] ]

  ndofs, el2Ms, eh1Ms, kopMs, iteMs, el2Ns, eh1Ns, kopNs, iteNs = conv_test(nes[1],1)

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

  ndofs, el2Ms, eh1Ms, kopMs, iteMs, el2Ns, eh1Ns, kopNs, iteNs = conv_test(nes[2],2)

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

  ndofs, el2Ms, eh1Ms, kopMs, iteMs, el2Ns, eh1Ns, kopNs, iteNs = conv_test(nes[3],3)

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

  # ndofs, el2Ms, eh1Ms, kopMs, iteMs, el2Ns, eh1Ns, kopNs, iteNs = conv_test(nes[4],4)

  # data = (log10.(ndofs),log10.(el2Ms))
  # lineML24 = lines!(axL,data,color=:red,linewidth=2)
  # scatML24 = scatter!(axL,data,color=:red,marker=:circle,markersize=6.0)

  # data = (log10.(ndofs),log10.(eh1Ms))
  # lineMH14 = lines!(axH,data,color=:red,linewidth=2)
  # scatMH14 = scatter!(axH,data,color=:red,marker=:xcross,markersize=6.0)

  # data = (log10.(ndofs),log10.(kopMs))
  # lineMCN4 = lines!(axC,data,color=:red,linewidth=2)
  # scatMCN4 = scatter!(axC,data,color=:red,marker=:diamond,markersize=6.0)

  # data = (log10.(ndofs),log10.(iteMs))
  # lineMIT4 = lines!(axI,data,color=:red,linewidth=2)
  # scatMIT4 = scatter!(axI,data,color=:red,marker=:utriangle,markersize=6.0)

  # # data = (log10.(ndofs),log10.(el2Ns))
  # # lineNL24 = lines!(axL,data,color=:blue,linestyle=:dash,linewidth=2)
  # # scatNL24 = scatter!(axL,data,color=:blue,marker=:circle,markersize=6.0)

  # # data = (log10.(ndofs),log10.(eh1Ns))
  # # lineNH14 = lines!(axH,data,color=:blue,linestyle=:dash,linewidth=2)
  # # scatNH14 = scatter!(axH,data,color=:blue,marker=:xcross,markersize=6.0)

  # data = (log10.(ndofs),log10.(kopNs))
  # lineNCN4 = lines!(axC,data,color=:blue,linestyle=:dash,linewidth=2)
  # scatNCN4 = scatter!(axC,data,color=:blue,marker=:diamond,markersize=6.0)

  # data = (log10.(ndofs),log10.(iteNs))
  # lineNIT4 = lines!(axI,data,color=:blue,linestyle=:dash,linewidth=2)
  # scatNIT4 = scatter!(axI,data,color=:blue,marker=:utriangle,markersize=6.0)

  # ndofs, el2Ms, eh1Ms, kopMs, iteMs, el2Ns, eh1Ns, kopNs, iteNs = conv_test(nes[5],5)

  # data = (log10.(ndofs),log10.(el2Ms))
  # lineML25 = lines!(axL,data,color=:red,linewidth=2)
  # scatML25 = scatter!(axL,data,color=:red,marker=:circle,markersize=6.0)

  # data = (log10.(ndofs),log10.(eh1Ms))
  # lineMH15 = lines!(axH,data,color=:red,linewidth=2)
  # scatMH15 = scatter!(axH,data,color=:red,marker=:xcross,markersize=6.0)

  # data = (log10.(ndofs),log10.(kopMs))
  # lineMCN5 = lines!(axC,data,color=:red,linewidth=2)
  # scatMCN5 = scatter!(axC,data,color=:red,marker=:diamond,markersize=6.0)

  # data = (log10.(ndofs),log10.(iteMs))
  # lineMIT5 = lines!(axI,data,color=:red,linewidth=2)
  # scatMIT5 = scatter!(axI,data,color=:red,marker=:utriangle,markersize=6.0)

  # # data = (log10.(ndofs),log10.(el2Ns))
  # # lineNL25 = lines!(axL,data,color=:blue,linestyle=:dash,linewidth=2)
  # # scatNL25 = scatter!(axL,data,color=:blue,marker=:circle,markersize=6.0)

  # # data = (log10.(ndofs),log10.(eh1Ns))
  # # lineNH15 = lines!(axH,data,color=:blue,linestyle=:dash,linewidth=2)
  # # scatNH15 = scatter!(axH,data,color=:blue,marker=:xcross,markersize=6.0)

  # data = (log10.(ndofs),log10.(kopNs))
  # lineNCN5 = lines!(axC,data,color=:blue,linestyle=:dash,linewidth=2)
  # scatNCN5 = scatter!(axC,data,color=:blue,marker=:diamond,markersize=6.0)

  # # data = (log10.(ndofs),log10.(iteNs))
  # # lineNIT5 = lines!(axI,data,color=:blue,linestyle=:dash,linewidth=2)
  # # scatNIT5 = scatter!(axI,data,color=:blue,marker=:utriangle,markersize=6.0)

  limits!(axL,0.65,2.85,-13.5,-1.5)
  axL.xticks = 0.75:0.5:2.75
  axL.yticks = -13.0:1.0:-2.0
  axL.xlabel="log10(N^(1/D))"
  axL.ylabel="log10(L2 Abs error)"

  limits!(axH,0.65,2.85,-9.5,-0.5)
  axH.xticks = 0.75:0.5:2.75
  axH.yticks = -9.0:1.0:-1.0
  axH.xlabel="log10(N^(1/D))"
  axH.ylabel="log10(H1 Abs error)"

  limits!(axC,0.65,2.85,-0.5,13.5)
  axC.xticks = 0.75:0.5:2.75
  axC.yticks = 0.0:1.0:13.0
  axC.xlabel="log10(N^(1/D))"
  axC.ylabel="log10(Condition number)"

  limits!(axI,0.65,2.85,-0.25,4.25)
  axI.xticks = 0.75:0.5:2.75
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

  # save("ModalC0AgFEM2D.png",scene)
  save("ModalC0AgFEM3D.png",scene)
  scene

end

plot()
