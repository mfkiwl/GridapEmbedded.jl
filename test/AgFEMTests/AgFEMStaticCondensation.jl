module AgFEMStaticCondensation

  using Gridap

  using MPI
  using GridapPETSc

  using GridapEmbedded
  using GridapEmbedded.AgFEM
  using Test

  using GridapEmbedded.Interfaces: CUT, IN
  using SparseArrays: SparseMatrixCSC

  include("StaticCondensation.jl")

  # using Profile
  # using BenchmarkTools

  if !MPI.Initialized()
    MPI.Init()
  end

  tol = 1e-16
  maxits = 1000
  GridapPETSc.Init(["-ksp_type", "cg",
                    "-ksp_rtol", "$tol",
                    "-ksp_max_it", "$maxits",
                    "-ksp_norm_type", "unpreconditioned",
                    "-pc_type","gamg",
                    "-pc_gamg_type","agg",
                    "-mg_coarse_sub_pc_type","cholesky",
                    "-mg_coarse_sub_pc_factor_mat_ordering_type","nd",
                    "-pc_gamg_process_eq_limit","50",
                    "-pc_gamg_square_graph","0",
                    "-pc_gamg_agg_nsmooths","1"])

  u(x) = (x[1]+x[2])^5
  g(x) = -Δ(u)(x)
  ud(x) = u(x)

  const R = 0.42

  n = 12
  order = 5

  geom = disk(R,x0=Point(0.5,0.5))
  partition = (n,n)
  domain = (0,1,0,1)

  bgmodel = CartesianDiscreteModel(domain,partition)
  h = (domain[2]-domain[1])/n
  cutgeo = cut(bgmodel,geom)
  model = DiscreteModel(cutgeo)

  strategy = AggregateAllCutCells()
  aggregates = aggregate(strategy,cutgeo,geom)
  bboxes = compute_cell_to_dface_bboxes(model,aggregates)

  Ωᵇ = Triangulation(bgmodel)
  Ω = Triangulation(cutgeo)
  Γ = EmbeddedBoundary(cutgeo)

  n_Γ = get_normal_vector(Γ)

  cutdeg, degree = 2*num_dims(model)*order, 2*order
  dΩ = Measure(Ω,cutdeg,degree)
  dΓ = Measure(Γ,cutdeg)

  reffe = ReferenceFE(modalC0,Float64,order,bboxes)
  Vstd = TestFESpace(model,reffe,conformity=:H1)
  V = AgFEMSpace(Vstd,aggregates,modalC0)
  U = TrialFESpace(V)

  γ = 5.0*order^2

  a(u,v) =
    ∫( ∇(v)⋅∇(u) )dΩ +
    ∫( (γ/h)*v*u  - v*(n_Γ⋅∇(u)) - (n_Γ⋅∇(v))*u )dΓ

  l(v) =
    ∫( v*g )dΩ +
    ∫( (γ/h)*v*ud - (n_Γ⋅∇(v))*ud )dΓ

  cell_matvec, bface_matvec = compute_contributions(U,V,a,l,Ω,Γ)
  cell_matvec, ccell_to_bgcell =
    combine_cell_and_bface_contribs(Ωᵇ,Ω,Γ,cell_matvec,bface_matvec)
  cell_matvec = attach_constraints_cols(U,cell_matvec,ccell_to_bgcell)
  cell_matvec = attach_constraints_rows(V,cell_matvec,ccell_to_bgcell)

  cell_dofs = lazy_map(Reindex(get_cell_dof_ids(V)),ccell_to_bgcell)
  idof_to_dof, dof_to_idof = compute_itfc_dofs(cell_dofs,V)
  cell_ildof_ldof, cell_cldof_ldof, cell_idofs =
    compute_itfc_cell_dof_arrays(cell_dofs,dof_to_idof)

  dict = Dict()
  cell_imatvec = lazy_map(StaticCondensationMap(dict),cell_matvec,cell_ildof_ldof,cell_cldof_ldof)

  S,f = assemble_schur_system(cell_imatvec,cell_idofs,idof_to_dof)
  solver = PETScSolver()
  x = solve(solver,S,f)
  y = apply_interior_correction(x,V,
                                idof_to_dof,dof_to_idof,cell_matvec,
                                cell_dofs,cell_ildof_ldof,cell_cldof_ldof,
                                dict)

  # @profile solve_with_static_condensation()
  # Profile.print(maxdepth=10)
  uh = FEFunction(V,y)
  # writevtk(Ω,"trian",cellfields=["uh"=>uh])

  e = u - uh

  l2(u) = sqrt(sum( ∫( u*u )dΩ ))
  h1(u) = sqrt(sum( ∫( u*u + ∇(u)⋅∇(u) )dΩ ))

  el2 = l2(e)
  eh1 = h1(e)
  ul2 = l2(uh)
  uh1 = h1(uh)

  @test el2/ul2 < 1.e-8
  @test eh1/uh1 < 1.e-7

  GridapPETSc.Finalize()

  # if MPI.Initialized() & !MPI.Finalized()
  #   MPI.Finalize()
  # end

end # module
