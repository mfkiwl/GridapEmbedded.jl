module AggregateBoundingBoxesTests

using Gridap
using Gridap.Geometry
using GridapEmbedded
using GridapEmbedded.AgFEM
using Test

u(x) = x[1] - x[2]
f(x) = -Δ(u)(x)
ud(x) = u(x)

const R = 0.42

geom = disk(R,x0=Point(0.5,0.5))
n = 5; partition = (n,n)
domain = (0,1,0,1)

# geom = sphere(R,x0=Point(0.5,0.5,0.5))
# n = 5; partition = (n,n,n)
# domain = (0,1,0,1,0,1)

bgmodel = CartesianDiscreteModel(domain,partition)
const h = (domain[2]-domain[1])/n
cutgeo = cut(bgmodel,geom)
model = DiscreteModel(cutgeo)

strategy = AggregateAllCutCells()
aggregates = aggregate(strategy,cutgeo,geom)
bboxes = compute_cell_to_dface_bboxes(model,aggregates)

Ω_bg = Triangulation(bgmodel)
Ω = Triangulation(cutgeo)
Γ = EmbeddedBoundary(cutgeo)

n_Γ = get_normal_vector(Γ)

order = 1
degree = 2*order*num_dims(model)
dΩ = Measure(Ω,degree)
dΓ = Measure(Γ,degree)

reffe = ReferenceFE(modalC0,Float64,order,bboxes)
Vstd = TestFESpace(model,reffe,conformity=:H1)
V = AgFEMSpace(Vstd,aggregates)
U = TrialFESpace(V)

const γd = 10.0

a(u,v) =
  ∫( ∇(v)⋅∇(u) ) * dΩ +
  ∫( (γd/h)*v*u  - v*(n_Γ⋅∇(u)) - (n_Γ⋅∇(v))*u ) * dΓ

l(v) =
  ∫( v*f ) * dΩ +
  ∫( (γd/h)*v*ud - (n_Γ⋅∇(v))*ud ) * dΓ

op = AffineFEOperator(a,l,U,V)
uh = solve(op)

e = u - uh

l2(u) = sqrt(sum( ∫( u*u )*dΩ ))
h1(u) = sqrt(sum( ∫( u*u + ∇(u)⋅∇(u) )*dΩ ))

el2 = l2(e)
eh1 = h1(e)
ul2 = l2(uh)
uh1 = h1(uh)

colors = color_aggregates(aggregates,bgmodel)
writevtk(Ω_bg,"trian",celldata=["cellin"=>aggregates,"color"=>colors])
writevtk(Ω_bg,"trian_S",nsubcells=10,cellfields=["uh"=>uh])
writevtk(Ω,"trian_O",cellfields=["uh"=>uh])

@test el2/ul2 < 1.e-8
@test eh1/uh1 < 1.e-7

# sktrian = BoundaryTriangulation(bgmodel,tags=collect(1:9))

# bbminsx = [ dbboxes[1][i][1].data[1] for i in 1:length(dbboxes[1]) ]
# bbminsy = [ dbboxes[1][i][1].data[2] for i in 1:length(dbboxes[1]) ]

# bbmaxsx = [ dbboxes[1][i][2].data[1] for i in 1:length(dbboxes[1]) ]
# bbmaxsy = [ dbboxes[1][i][2].data[2] for i in 1:length(dbboxes[1]) ]

# writevtk(sktrian,"sktrian",
#          celldata=["aminsx"=>bbminsx,"aminsy"=>bbminsy,
#                    "bmaxsx"=>bbmaxsx,"bmaxsy"=>bbmaxsy])

# sktrian = BoundaryTriangulation(bgmodel,tags=collect(1:27))

# bbminsx = [ dbboxes[1][i][1].data[1] for i in 1:length(dbboxes[2]) ]
# bbminsy = [ dbboxes[1][i][1].data[2] for i in 1:length(dbboxes[2]) ]
# bbminsz = [ dbboxes[1][i][1].data[3] for i in 1:length(dbboxes[2]) ]

# bbmaxsx = [ dbboxes[1][i][2].data[1] for i in 1:length(dbboxes[2]) ]
# bbmaxsy = [ dbboxes[1][i][2].data[2] for i in 1:length(dbboxes[2]) ]
# bbmaxsz = [ dbboxes[1][i][2].data[3] for i in 1:length(dbboxes[2]) ]

# writevtk(sktrian,"sktrian",
#          celldata=["aminsx"=>bbminsx,"aminsy"=>bbminsy,"aminsz"=>bbminsz,
#                    "bmaxsx"=>bbmaxsx,"bmaxsy"=>bbmaxsy,"bmaxsz"=>bbmaxsz])

end # module
