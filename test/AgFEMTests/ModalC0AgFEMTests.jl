module ModalC0AgFEMTests

  using DrWatson
  @quickactivate "GridapEmbedded"

  using Gridap
  using GridapEmbedded
  using GridapEmbedded.AgFEM

  using LinearAlgebra: cond
  using SparseArrays: SparseMatrixCSC

  using GridapPETSc
  using MPI

  const tol = 1e-16
  const maxits = 10000
  const R = 0.42
  const ω = 2*pi
  exact(x) = x[1]+x[2]
  sinus(x) = sin(ω*x[1])*x[2]

  function compute(n::Int,k::Int,d::Int,t::Int,s::Int,g::Int)

    u(x) = (1-s)*exact(x)+s*sinus(x)
    f(x) = -Δ(u)(x)
    ud(x) = u(x)

    if d == 2
      if g == 0
        geom = disk(R,x0=Point(0.5,0.5))
      else
        geom = square(L=0.63,x0=Point(0.5,0.5))
      end
      partition = (n,n)
      domain = (0,1,0,1)
    else
      if g == 0
        geom = sphere(R,x0=Point(0.5,0.5,0.5))
      else
        geom = doughnut(0.21,0.1,x0=Point(0.5,0.5,0.5))
      end
      partition = (n,n,n)
      domain = (0,1,0,1,0,1)
    end

    bgmodel = CartesianDiscreteModel(domain,partition)
    h = (domain[2]-domain[1])/n
    cutgeo = cut(bgmodel,geom)
    model = DiscreteModel(cutgeo)

    if t == 0
      γ₀ = 0.65
      eᵧ = (3-2/num_dims(model))/(2*k+1-2/num_dims(model))
      γ = γ₀^eᵧ
      strategy = AggregateCutCellsByThreshold(γ)
    else
      strategy = AggregateAllCutCells()
    end

    aggregates = aggregate(strategy,cutgeo,geom)
    bboxes = compute_cell_to_dface_bboxes(model,aggregates)

    Ω_bg = Triangulation(bgmodel)
    Ω = Triangulation(cutgeo)
    Γ = EmbeddedBoundary(cutgeo)

    n_Γ = get_normal_vector(Γ)

    D = num_dims(model)
    if k != 1 
      cdegm, cdegs, deg = 2*D*k, 2*D*(k-1), 2*k
    else
      cdegm, cdegs, deg = 2*D*k, 2*D*k, 2*k
    end
    dΩ = Measure(Ω,cdegs,deg)
    dO = Measure(Ω,cdegm,deg)
    dΓ = Measure(Γ,cdegm)

    γd = 5.0*k^2

    a(u,v) =
      ∫( ∇(u)⋅∇(v) ) * dΩ +
      ∫( (γd/h)*u*v  - v*(n_Γ⋅∇(u)) - (n_Γ⋅∇(v))*u ) * dΓ

    l(v) =
      ∫( v*f ) * dΩ +
      ∫( (γd/h)*v*ud - (n_Γ⋅∇(v))*ud ) * dΓ

    l2(u) = sqrt(sum( ∫( u*u )*dO ))
    h1(u) = sqrt(sum( ∫( ∇(u)⋅∇(u) )*dΩ ))

    reffe = ReferenceFE(modalC0,Float64,k,bboxes)
    Vstd = TestFESpace(model,reffe,conformity=:H1)
    V = AgFEMSpace(Vstd,aggregates,modalC0)
    U = TrialFESpace(V)

    sop = @timed op = AffineFEOperator(a,l,U,V)
    scn = @timed kopM = cond(get_matrix(op),1)

    ls = PETScSolver()
    solver = LinearFESolver(ls)

    sso = @timed uh = solve(solver,op)
    iteM = PETSc_get_number_of_iterations(ls)

    e = u - uh
    el2M = l2(e)
    eh1M = h1(e)

    if k != 1
      reffe = ReferenceFE(lagrangian,Float64,k)
      Vstd = TestFESpace(model,reffe,conformity=:H1)
      V = AgFEMSpace(Vstd,aggregates,lagrangian)
      U = TrialFESpace(V)

      op = AffineFEOperator(a,l,U,V)
      kopN = cond(get_matrix(op),1)

      ls = PETScSolver()
      solver = LinearFESolver(ls)

      uh = solve(solver,op)
      iteN = PETSc_get_number_of_iterations(ls)

      e = u - uh
      el2N = l2(e)
      eh1N = h1(e)
    else
      el2N = el2M
      eh1N = eh1M
      kopN = kopM
      iteN = iteM
    end

    num_free_dofs(U)^(1/D), el2M, eh1M, kopM, iteM, el2N, eh1N, kopN, iteN, sop.time, scn.time, sso.time, sop.bytes, scn.bytes, sso.bytes

  end

  function compute(case::Dict)
    @unpack n,k,d,t,s,g = case
    udofs, el2M, eh1M, kopM, iteM, el2N, eh1N, kopN, iteN, top, tcn, tso, bop, bcn, bso = compute(n,k,d,t,s,g)
    results = @dict udofs el2M eh1M kopM iteM el2N eh1N kopN iteN top tcn tso bop bcn bso
    merge(case,results)
  end

  function compute_and_save(case::Dict)
    produce_or_load(
      datadir(),
      case,
      compute,
      prefix="modalC0",
      tag=true,
      verbose=true
    )
    return true
  end

  function compute()

    if !MPI.Initialized()
      MPI.Init()
    end

    GridapPETSc.Init(["-ksp_type","cg",
                      "-ksp_rtol","$tol",
                      "-ksp_max_it","$maxits",
                      "-ksp_norm_type","unpreconditioned",
                      "-pc_type","gamg",
                      "-pc_gamg_type","agg",
                      "-pc_gamg_esteig_ksp_type","cg",
                      "-mg_levels_esteig_ksp_type","cg",
                      "-mg_coarse_sub_pc_type","cholesky",
                      "-mg_coarse_sub_pc_factor_mat_ordering_type","nd",
                      "-pc_gamg_process_eq_limit","50",
                      "-pc_gamg_square_graph","0",
                      "-pc_gamg_agg_nsmooths","1"])

    @info "Training"
    compute(6,3,2,1,0,0)
    @info "Producing"
    params = Dict(
      :n => [6,8,10,12,14,16,18,20,22],
      :k => [3],
      :d => [2],
      :t => [1],
      :s => [0],
      :g => [0]
    )
    dicts = dict_list(params)
    map(compute_and_save,dicts)

    GridapPETSc.Finalize()

    if MPI.Initialized() & !MPI.Finalized()
      MPI.Finalize()
    end

  end

  # compute()
  export tol, maxits
  export compute, compute_and_save

end
