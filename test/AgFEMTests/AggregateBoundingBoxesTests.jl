module AggregateBoundingBoxesTests

  using Gridap
  using GridapEmbedded
  using GridapEmbedded.AgFEM
  # using GridapPETSc
  # using MPI
  using Test

  # if !MPI.Initialized()
  #   MPI.Init()
  # end

  # tol = 1e-16
  # maxits = 1000
  # GridapPETSc.Init(["-ksp_type", "cg",
  #                   "-ksp_monitor",
  #                   "-ksp_rtol", "$tol",
  #                   "-ksp_converged_reason",
  #                   "-ksp_max_it", "$maxits",
  #                   "-ksp_norm_type", "unpreconditioned",
  #                   "-ksp_view",
  #                   "-pc_type","gamg",
  #                   "-pc_gamg_type","agg",
  #                   "-pc_gamg_esteig_ksp_type","cg",
  #                   "-mg_levels_esteig_ksp_type","cg",
  #                   "-mg_coarse_sub_pc_type","cholesky",
  #                   "-mg_coarse_sub_pc_factor_mat_ordering_type","nd",
  #                   "-pc_gamg_process_eq_limit","50",
  #                   "-pc_gamg_square_graph","0",
  #                   "-pc_gamg_agg_nsmooths","1"])

  # using AlgebraicMultigrid
  # import IterativeSolvers: cg

  u(x) = x[1]+x[2]
  f(x) = -Δ(u)(x)
  ud(x) = u(x)

  const R = 0.42

  n = 12
  order = 1

  # geom = square(L=0.63,x0=Point(0.5,0.5))
  geom = disk(R,x0=Point(0.5,0.5))
  partition = (n,n)
  domain = (0,1,0,1)

  # geom = sphere(R,x0=Point(0.5,0.5,0.5))
  # geom = doughnut(0.21,0.10,x0=Point(0.5,0.5,0.5))
  # partition = (n,n,n)
  # domain = (0,1,0,1,0,1)

  bgmodel = CartesianDiscreteModel(domain,partition)
  h = (domain[2]-domain[1])/n
  cutgeo = cut(bgmodel,geom)
  model = DiscreteModel(cutgeo)

  γ₀ = 0.65
  eᵧ = (3-2/num_dims(model))/(2*order+1-2/num_dims(model))
  γ = γ₀^eᵧ
  strategy = AggregateCutCellsByThreshold(γ)
  # strategy = AggregateAllCutCells()
  aggregates = aggregate(strategy,cutgeo,geom)
  bboxes = compute_cell_to_dface_bboxes(model,aggregates)

  Ω_bg = Triangulation(bgmodel)
  Ω = Triangulation(cutgeo)
  Γ = EmbeddedBoundary(cutgeo)

  n_Γ = get_normal_vector(Γ)

  cutdeg, degree = 2*num_dims(model)*order, 2*order
  dΩ = Measure(Ω,cutdeg,degree)
  dΓ = Measure(Γ,cutdeg)

  # reffe = ReferenceFE(lagrangian,Float64,order)
  reffe = ReferenceFE(modalC0,Float64,order,bboxes)
  Vstd = TestFESpace(model,reffe,conformity=:H1)
  # V = AgFEMSpace(Vstd,aggregates,lagrangian)
  V = AgFEMSpace(Vstd,aggregates,modalC0)
  U = TrialFESpace(V)

  γd = 5.0*order^2

  a(u,v) =
    ∫( ∇(v)⋅∇(u) ) * dΩ +
    ∫( (γd/h)*v*u  - v*(n_Γ⋅∇(u)) - (n_Γ⋅∇(v))*u ) * dΓ

  l(v) =
    ∫( v*f ) * dΩ +
    ∫( (γd/h)*v*ud - (n_Γ⋅∇(v))*ud ) * dΓ

  op = AffineFEOperator(a,l,U,V)

  uh = solve(op)

  # A = get_matrix(op)
  # b = get_vector(op)
  # p = aspreconditioner(smoothed_aggregation(A))
  # ur, ch = cg(A,b,Pl=p,log=true)
  # uh = FEFunction(U,ur)

  # ass = SparseMatrixAssembler(SparseMatrixCSR{0,PetscReal,PetscInt},U,V)
  # op = AffineFEOperator(a,l,ass)

  # ls = PETScSolver()
  # solver = LinearFESolver(ls)

  # uh = solve(solver,op)
  # iters = PETSc_get_number_of_iterations(ls)

  e = u - uh

  l2(u) = sqrt(sum( ∫( u*u )*dΩ ))
  h1(u) = sqrt(sum( ∫( u*u + ∇(u)⋅∇(u) )*dΩ ))

  el2 = l2(e)
  eh1 = h1(e)
  ul2 = l2(uh)
  uh1 = h1(uh)

  # colors = color_aggregates(aggregates,bgmodel)
  # writevtk(Ω_bg,"trian",celldata=["cellin"=>aggregates,"color"=>colors])
  # writevtk(Ω_bg,"trian_S",nsubcells=10,cellfields=["uh"=>uh])
  # writevtk(Ω,"trian_O",cellfields=["uh"=>uh])

  # GridapPETSc.Finalize()

  @test el2/ul2 < 1.e-8
  @test eh1/uh1 < 1.e-7

  # if MPI.Initialized() & !MPI.Finalized()
  #   MPI.Finalize()
  # end

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
