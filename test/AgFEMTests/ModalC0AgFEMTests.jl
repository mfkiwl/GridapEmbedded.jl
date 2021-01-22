using Gridap
using Gridap.Fields
using Gridap.Geometry
using GridapEmbedded
using GridapEmbedded.AgFEM
using LinearAlgebra: cond

using Makie
using AbstractPlotting
using AbstractPlotting.MakieLayout

u(x) = (x[1] - x[2])^5
f(x) = -Δ(u)(x)
ud(x) = u(x)

const R = 0.4

function run_modal(n,order)

  geom = disk(R,x0=Fields.Point(0.5,0.5))
  partition = (n,n)
  domain = (0,1,0,1)

  # geom = sphere(R,x0=Fields.Point(0.0,0.0,0.0))
  # n = 5; partition = (n,n,n)
  # domain = (0,1,0,1,0,1)

  bgmodel = CartesianDiscreteModel(domain,partition)
  h = (domain[2]-domain[1])/n
  cutgeo = cut(bgmodel,geom)
  model = DiscreteModel(cutgeo)

  strategy = AggregateAllCutCells()
  aggregates = aggregate(strategy,cutgeo,geom)
  bboxes = compute_cell_to_dface_bboxes(model,aggregates)

  Ω_bg = Triangulation(bgmodel)
  Ω = Triangulation(cutgeo)
  Γ = EmbeddedBoundary(cutgeo)

  n_Γ = get_normal_vector(Γ)

  D = num_dims(model)
  cutdeg, degree = 2*order*D, 2*order
  dΩ = Measure(Ω,cutdeg,degree)
  dΓ = Measure(Γ,cutdeg)

  reffe = ReferenceFE(modalC0,Float64,order,bboxes)
  Vstd = TestFESpace(model,reffe,conformity=:H1)
  V = AgFEMSpace(Vstd,aggregates,modalC0)
  U = TrialFESpace(V)

  γd = 2.5*order^2

  a(u,v) =
    ∫( ∇(v)⋅∇(u) ) * dΩ +
    ∫( (γd/h)*v*u  - v*(n_Γ⋅∇(u)) - (n_Γ⋅∇(v))*u ) * dΓ

  l(v) =
    ∫( v*f ) * dΩ +
    ∫( (γd/h)*v*ud - (n_Γ⋅∇(v))*ud ) * dΓ

  op = AffineFEOperator(a,l,U,V)
  kop = cond(Array(get_matrix(op)))
  uh = solve(op)

  e = u - uh

  l2(u) = sqrt(sum( ∫( u*u )*dΩ ))
  h1(u) = sqrt(sum( ∫( u*u + ∇(u)⋅∇(u) )*dΩ ))

  el2 = l2(e)
  eh1 = h1(e)

  num_free_dofs(U)^(1/D), el2, eh1, kop

end

function run_nodal(n,order)

  geom = disk(R,x0=Fields.Point(0.5,0.5))
  partition = (n,n)
  domain = (0,1,0,1)

  # geom = sphere(R,x0=Fields.Point(0.0,0.0,0.0))
  # n = 5; partition = (n,n,n)
  # domain = (0,1,0,1,0,1)

  bgmodel = CartesianDiscreteModel(domain,partition)
  h = (domain[2]-domain[1])/n
  cutgeo = cut(bgmodel,geom)
  model = DiscreteModel(cutgeo)

  strategy = AggregateAllCutCells()
  aggregates = aggregate(strategy,cutgeo,geom)

  Ω_bg = Triangulation(bgmodel)
  Ω = Triangulation(cutgeo)
  Γ = EmbeddedBoundary(cutgeo)

  n_Γ = get_normal_vector(Γ)

  D = num_dims(model)
  cutdeg, degree = 2*order*D, 2*order
  dΩ = Measure(Ω,cutdeg,degree)
  dΓ = Measure(Γ,cutdeg)

  reffe = ReferenceFE(lagrangian,Float64,order)
  Vstd = TestFESpace(model,reffe,conformity=:H1)
  V = AgFEMSpace(Vstd,aggregates,lagrangian)
  U = TrialFESpace(V)

  γd = 2.5*order^2

  a(u,v) =
    ∫( ∇(v)⋅∇(u) ) * dΩ +
    ∫( (γd/h)*v*u  - v*(n_Γ⋅∇(u)) - (n_Γ⋅∇(v))*u ) * dΓ

  l(v) =
    ∫( v*f ) * dΩ +
    ∫( (γd/h)*v*ud - (n_Γ⋅∇(v))*ud ) * dΓ

  op = AffineFEOperator(a,l,U,V)
  kop = cond(Array(get_matrix(op)))
  uh = solve(op)

  e = u - uh

  l2(u) = sqrt(sum( ∫( u*u )*dΩ ))
  h1(u) = sqrt(sum( ∫( u*u + ∇(u)⋅∇(u) )*dΩ ))

  el2 = l2(e)
  eh1 = h1(e)

  num_free_dofs(U)^(1/D), el2, eh1, kop

end

function conv_test(ns,order,run)

  ndfs = Float64[]
  el2s = Float64[]
  eh1s = Float64[]
  kops = Float64[]

  for n in ns

    ndf, el2, eh1, kop = run(n,order)

    push!(ndfs,ndf)
    push!(el2s,el2)
    push!(eh1s,eh1)
    push!(kops,kop)

  end

  (ndfs, el2s, eh1s, kops)

end

outer_padding = 30
scene, layout = layoutscene(outer_padding, resolution = (1000, 400),
                            backgroundcolor = RGBf0(0.99, 0.99, 0.99))

axL = layout[1,1] = LAxis(scene)
axH = layout[1,2] = LAxis(scene)
axC = layout[1,3] = LAxis(scene)

nes = [ [8,16,32,64,128], [6,12,24,48], [6,12,24,36], [6,12,18,24] ]
# nes = [ [8,32], [6,48], [6,36], [6,24] ]

ndfs, el2s, eh1s, kops = conv_test(nes[1],1,run_modal)

data = (log10.(ndfs),log10.(el2s))
lineML21 = lines!(axL,data,color=:red,linewidth=2)
scatML21 = scatter!(axL,data,color=:red,marker=:circle,markersize=6.0)

data = (log10.(ndfs),log10.(eh1s))
lineMH11 = lines!(axH,data,color=:red,linewidth=2)
scatMH11 = scatter!(axH,data,color=:red,marker=:xcross,markersize=6.0)

data = (log10.(ndfs),log10.(kops))
lineMCN1 = lines!(axC,data,color=:red,linewidth=2)
scatMCN1 = scatter!(axC,data,color=:red,marker=:diamond,markersize=6.0)

ndfs, el2s, eh1s, kops = conv_test(nes[1],1,run_nodal)

data = (log10.(ndfs),log10.(el2s))
lineNL21 = lines!(axL,data,color=:blue,linestyle=:dash,linewidth=2)
scatNL21 = scatter!(axL,data,color=:blue,marker=:circle,markersize=6.0)

data = (log10.(ndfs),log10.(eh1s))
lineNH11 = lines!(axH,data,color=:blue,linestyle=:dash,linewidth=2)
scatNH11 = scatter!(axH,data,color=:blue,marker=:xcross,markersize=6.0)

data = (log10.(ndfs),log10.(kops))
lineNCN1 = lines!(axC,data,color=:blue,linestyle=:dash,linewidth=2)
scatNCN1 = scatter!(axC,data,color=:blue,marker=:diamond,markersize=6.0)

ndfs, el2s, eh1s, kops = conv_test(nes[2],2,run_modal)

data = (log10.(ndfs),log10.(el2s))
lineML22 = lines!(axL,data,color=:red,linewidth=2)
scatML22 = scatter!(axL,data,color=:red,marker=:circle,markersize=6.0)

data = (log10.(ndfs),log10.(eh1s))
lineMH12 = lines!(axH,data,color=:red,linewidth=2)
scatMH12 = scatter!(axH,data,color=:red,marker=:xcross,markersize=6.0)

data = (log10.(ndfs),log10.(kops))
lineMCN2 = lines!(axC,data,color=:red,linewidth=2)
scatMCN2 = scatter!(axC,data,color=:red,marker=:diamond,markersize=6.0)

ndfs, el2s, eh1s, kops = conv_test(nes[2],2,run_nodal)

data = (log10.(ndfs),log10.(el2s))
lineNL22 = lines!(axL,data,color=:blue,linestyle=:dash,linewidth=2)
scatNL22 = scatter!(axL,data,color=:blue,marker=:circle,markersize=6.0)

data = (log10.(ndfs),log10.(eh1s))
lineNH12 = lines!(axH,data,color=:blue,linestyle=:dash,linewidth=2)
scatNH12 = scatter!(axH,data,color=:blue,marker=:xcross,markersize=6.0)

data = (log10.(ndfs),log10.(kops))
lineNCN2 = lines!(axC,data,color=:blue,linestyle=:dash,linewidth=2)
scatNCN2 = scatter!(axC,data,color=:blue,marker=:diamond,markersize=6.0)

ndfs, el2s, eh1s, kops = conv_test(nes[3],3,run_modal)

data = (log10.(ndfs),log10.(el2s))
lineML23 = lines!(axL,data,color=:red,linewidth=2)
scatML23 = scatter!(axL,data,color=:red,marker=:circle,markersize=6.0)

data = (log10.(ndfs),log10.(eh1s))
lineMH13 = lines!(axH,data,color=:red,linewidth=2)
scatMH13 = scatter!(axH,data,color=:red,marker=:xcross,markersize=6.0)

data = (log10.(ndfs),log10.(kops))
lineMCN3 = lines!(axC,data,color=:red,linewidth=2)
scatMCN3 = scatter!(axC,data,color=:red,marker=:diamond,markersize=6.0)

ndfs, el2s, eh1s, kops = conv_test(nes[3],3,run_nodal)

data = (log10.(ndfs),log10.(el2s))
lineNL23 = lines!(axL,data,color=:blue,linestyle=:dash,linewidth=2)
scatNL23 = scatter!(axL,data,color=:blue,marker=:circle,markersize=6.0)

data = (log10.(ndfs),log10.(eh1s))
lineNH13 = lines!(axH,data,color=:blue,linestyle=:dash,linewidth=2)
scatNH13 = scatter!(axH,data,color=:blue,marker=:xcross,markersize=6.0)

data = (log10.(ndfs),log10.(kops))
lineNCN3 = lines!(axC,data,color=:blue,linestyle=:dash,linewidth=2)
scatNCN3 = scatter!(axC,data,color=:blue,marker=:diamond,markersize=6.0)

ndfs, el2s, eh1s, kops = conv_test(nes[4],4,run_modal)

data = (log10.(ndfs),log10.(el2s))
lineML24 = lines!(axL,data,color=:red,linewidth=2)
scatML24 = scatter!(axL,data,color=:red,marker=:circle,markersize=6.0)

data = (log10.(ndfs),log10.(eh1s))
lineMH14 = lines!(axH,data,color=:red,linewidth=2)
scatMH14 = scatter!(axH,data,color=:red,marker=:xcross,markersize=6.0)

data = (log10.(ndfs),log10.(kops))
lineMCN4 = lines!(axC,data,color=:red,linewidth=2)
scatMCN4 = scatter!(axC,data,color=:red,marker=:diamond,markersize=6.0)

ndfs, el2s, eh1s, kops = conv_test(nes[4],4,run_nodal)

data = (log10.(ndfs),log10.(el2s))
lineNL24 = lines!(axL,data,color=:blue,linestyle=:dash,linewidth=2)
scatNL24 = scatter!(axL,data,color=:blue,marker=:circle,markersize=6.0)

data = (log10.(ndfs),log10.(eh1s))
lineNH14 = lines!(axH,data,color=:blue,linestyle=:dash,linewidth=2)
scatNH14 = scatter!(axH,data,color=:blue,marker=:xcross,markersize=6.0)

data = (log10.(ndfs),log10.(kops))
lineNCN4 = lines!(axC,data,color=:blue,linestyle=:dash,linewidth=2)
scatNCN4 = scatter!(axC,data,color=:blue,marker=:diamond,markersize=6.0)

limits!(axL,0.65,2.1,-9.5,-1.5)
axL.xticks = 0.75:0.25:2.0
axL.yticks = -9.0:1.0:-2.0
axL.xlabel="log10(N)"
axL.ylabel="log10(Abs error)"

limits!(axH,0.65,2.1,-7.5,-0.5)
axH.xticks = 0.75:0.25:2.0
axH.yticks = -7.0:1.0:-1.0
axH.xlabel="log10(N)"
axH.ylabel="log10(Abs error)"

limits!(axC,0.65,2.1,-0.5,9.5)
axC.xticks = 0.75:0.25:2.0
axC.yticks = 0.0:1.0:9.0
axC.xlabel="log10(N)"
axC.ylabel="log10(Condition number)"

axL.xticksize = 2.0; axL.yticksize = 2.0
axL.xticklabelsize = 12.0; axL.yticklabelsize = 12.0

axH.xticksize = 2.0; axH.yticksize = 2.0
axH.xticklabelsize = 12.0; axH.yticklabelsize = 12.0

axC.xticksize = 2.0; axC.yticksize = 2.0
axC.xticklabelsize = 12.0; axC.yticklabelsize = 12.0

mark1 = MarkerElement(color=:gray,marker=:circle,strokecolor=:black)
mark2 = MarkerElement(color=:gray,marker=:xcross,strokecolor=:black)
mark3 = MarkerElement(color=:gray,marker=:diamond,strokecolor=:black)
legmarkers = [ lineML21, lineNL21, mark1, mark2, mark3 ]
legnames = [ "Modal", "Nodal", "L2(e)", "H1(e)", "κ(A)" ]
layout[1,4] = LLegend( scene, legmarkers, legnames, orientation = :vertical )

save("ModalC0AgFEM2D.png",scene)
scene
