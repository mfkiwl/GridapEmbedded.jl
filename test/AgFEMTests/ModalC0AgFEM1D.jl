module ModalC0AgFEM1D

using Gridap
using Gridap.Arrays
using Gridap.Geometry
using GridapEmbedded
using GridapEmbedded.AgFEM
using Gridap.FESpaces
using Gridap.ReferenceFEs
using LinearAlgebra: cond
using Test

const maxod = 2
u(x) = -x[1]^maxod+maxod*x[1]
g(x) = ∇(u)(x)
f(x) = -Δ(u)(x)

hs(x::Real) = x > 1.0 ? 1.0 : 0.0

function stretching(x::Point)
  m = zeros(length(x))
  m[1] = x[1] + hs(x[1])*κ
  Point(m)
end

κ = 0.5
n = 2; partition = (n); domain = (0,2)
model = CartesianDiscreteModel(domain,partition,map=stretching)

aggregates = Int32[1,1]
bboxes = [Point{1,Float64}[(0.0),(1.0),(0.0),(1.0),(0.0),(2.0+κ)],
          Point{1,Float64}[(0.0),(1.0),(0.0),(1.0),(0.0),(1.0)]]

Ω = Triangulation(model)
Γ = BoundaryTriangulation(model)
n_Γ = get_normal_vector(Γ)
# writevtk(Ω,"debug")
od = 3; degree = 2*od
dΩ = Measure(Ω,degree)
dΓ = Measure(Γ,degree)

# reffe = ReferenceFE(lagrangian,Float64,od)
reffe = ReferenceFE(modalC0,Float64,od,bboxes)
Vstd = TestFESpace(model,reffe,dirichlet_tags=[1])

aggdof_to_dof = vcat(2,collect(2+od:2*od))
aggdof_to_dofs_ptrs = [1+i*(od+1) for i in 0:od]
aggdof_to_dofs_data = repeat(vcat(-1,1,collect(3:od+1)),od)
aggdof_to_dofs = Table(aggdof_to_dofs_data,aggdof_to_dofs_ptrs)
sh1 = get_cell_shapefuns(Vstd).cell_basis[1]
df2 = get_cell_dof_basis(Vstd).cell_dof[2]
ld2 = ( ( 1.0 + κ ) .* df2.dof_basis.nodes ) .+ 1.0
ld2 = LagrangianDofBasis(df2.dof_basis,ld2)
df2to1 = linear_combination(df2,ld2)
aggdof_to_coeffs_data = transpose(evaluate(df2to1,sh1))[od+2:end]
aggdof_to_coeffs = Table(aggdof_to_coeffs_data,aggdof_to_dofs_ptrs)

V = FESpaceWithLinearConstraints(aggdof_to_dof,aggdof_to_dofs,aggdof_to_coeffs,Vstd)
# V = Vstd
U = TrialFESpace(V,u)

a(u,v) = ∫( ∇(u)⋅∇(v) )*dΩ
l(v) = ∫( f*v )*dΩ + ∫( (g⋅n_Γ)*v )*dΓ

op = AffineFEOperator(a,l,U,V)
kop = cond(Array(get_matrix(op)))
uh = solve(op)

e = u - uh

l2(u) = sqrt(sum( ∫( u*u )*dΩ ))
h1(u) = sqrt(sum( ∫( u*u + ∇(u)⋅∇(u) )*dΩ ))

el2 = l2(e)
eh1 = h1(e)
ul2 = l2(uh)
uh1 = h1(uh)

# writevtk(Ω,"trian",nsubcells=10,cellfields=["uh"=>uh])

# @test el2/ul2 < 1.e-8
# @test eh1/uh1 < 1.e-7

end # module
