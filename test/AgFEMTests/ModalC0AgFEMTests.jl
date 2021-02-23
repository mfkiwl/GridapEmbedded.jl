using Gridap
using Gridap.Fields
using Gridap.Geometry
using GridapEmbedded
using GridapEmbedded.AgFEM
using LinearAlgebra: cond
using SparseArrays: SparseMatrixCSC

using GridapPETSc
using MPI

if !MPI.Initialized()
  MPI.Init()
end

tol = 1e-16
maxits = 10000
GridapPETSc.Init(["-ksp_type", "cg",
                  "-ksp_rtol", "$tol",
                  "-ksp_max_it", "$maxits",
                  "-ksp_norm_type", "unpreconditioned",
                  "-pc_type","gamg",
                  "-pc_gamg_type","agg",
                  "-pc_gamg_esteig_ksp_type","cg",
                  "-mg_levels_esteig_ksp_type","cg",
                  "-mg_coarse_sub_pc_type","cholesky",
                  "-mg_coarse_sub_pc_factor_mat_ordering_type","nd",
                  "-pc_gamg_process_eq_limit","50",
                  "-pc_gamg_square_graph","0",
                  "-pc_gamg_agg_nsmooths","1"])

const R = 0.42

# u(x) = x[1] + x[2]
const k = 2*pi
u(x) = sin(k*x[1])*x[2]
f(x) = -Δ(u)(x)
ud(x) = u(x)

function run(n::Int,order::Int)

  geom = disk(R,x0=Fields.Point(0.5,0.5))
  partition = (n,n)
  domain = (0,1,0,1)

  # geom = sphere(R,x0=Fields.Point(0.5,0.5,0.5))
  # partition = (n,n,n)
  # domain = (0,1,0,1,0,1)

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
  cutdeg, degree = 2*D*order, 2*order
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

  ass = SparseMatrixAssembler(SparseMatrixCSR{0,PetscReal,PetscInt},U,V)
  op = AffineFEOperator(a,l,ass)
  A = convert(SparseMatrixCSC, get_matrix(op))
  kopM = cond(A,1) # https://github.com/JuliaLang/julia/issues/6485

  ls = PETScSolver()
  solver = LinearFESolver(ls)

  uh = solve(solver,op)
  iteM = PETSc_get_number_of_iterations(ls)

  e = u - uh
  el2M = l2(e)
  eh1M = h1(e)

  reffe = ReferenceFE(lagrangian,Float64,order)
  Vstd = TestFESpace(model,reffe,conformity=:H1)
  V = AgFEMSpace(Vstd,aggregates,lagrangian)
  U = TrialFESpace(V)

  ass = SparseMatrixAssembler(SparseMatrixCSR{0,PetscReal,PetscInt},U,V)
  A = convert(SparseMatrixCSC, get_matrix(op))
  kopN = cond(A,1) # https://github.com/JuliaLang/julia/issues/6485

  ls = PETScSolver()
  solver = LinearFESolver(ls)

  uh = solve(solver,op)
  iteN = PETSc_get_number_of_iterations(ls)

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

  # nes = [ [8,16,32,64,128,256,512], [6,12,24,48,96,192], [6,12,24,48,96], [6,12,24,48,96], [6,12,24,48] ]
  nes = [ [8,16], [6,12,24], [6], [6,12,24], [6,12,24,48] ]

  ndofs, el2Ms, eh1Ms, kopMs, iteMs, el2Ns, eh1Ns, kopNs, iteNs = conv_test(nes[1],1)

end

plot()

GridapPETSc.Finalize()

if MPI.Initialized() & !MPI.Finalized()
  MPI.Finalize()
end
