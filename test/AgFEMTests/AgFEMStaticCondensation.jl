module AgFEMStaticCondensation

  using Gridap
  using Gridap.Geometry
  using Gridap.CellData
  using Gridap.Arrays
  using Gridap.Algebra
  using Gridap.FESpaces
  using SparseArrays

  using GridapEmbedded
  using GridapEmbedded.AgFEM
  using Test

  using GridapEmbedded.Interfaces: CUT, IN

  u(x) = (x[1]+x[2])^4
  g(x) = -Δ(u)(x)
  ud(x) = u(x)

  const R = 0.42

  n = 12
  order = 4

  geom = disk(R,x0=Point(0.5,0.5))
  partition = (n,n)
  domain = (0,1,0,1)

  bgmodel = CartesianDiscreteModel(domain,partition)
  h = (domain[2]-domain[1])/n
  cutgeo = cut(bgmodel,geom)
  model = DiscreteModel(cutgeo)

  active_model = DiscreteModel(cutgeo,geom,(CUT,IN))
  Ωᵃ = Triangulation(active_model)

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

  # Build cell and bface contributions
  du = get_cell_shapefuns_trial(V)
  dv = get_cell_shapefuns(V)

  cell_mat = a(du,dv)[Ω]
  cell_vec = l(dv)[Ω]
  cell_matvec = pair_arrays(cell_mat,cell_vec)

  bface_mat = a(du,dv)[Γ]
  bface_vec = l(dv)[Γ]

  # Compress cell and face contributions
  cell_to_bgcell = get_cell_to_bgcell(Ω)
  ccell_to_first_cell = compress(cell_to_bgcell)
  ccell_to_matvec, ccell_to_bgcell = compress(cell_matvec,cell_to_bgcell,ccell_to_first_cell)

  function compress_to_bgcells(t::Triangulation,bg::Triangulation)
    cell_to_bgcell = get_cell_to_bgcell(t)
    bgcell_to_cells_ptrs = zeros(Int32,num_cells(bg)+1)
    for bgcell in cell_to_bgcell
      bgcell_to_cells_ptrs[bgcell+1] += 1
    end
    length_to_ptrs!(bgcell_to_cells_ptrs)
    ndata = bgcell_to_cells_ptrs[end]-1
    bgcell_to_cells_data = zeros(Int32,ndata)
    for (cell,bgcell) in enumerate(cell_to_bgcell)
      p = bgcell_to_cells_ptrs[bgcell]
      bgcell_to_cells_data[p] = cell
      bgcell_to_cells_ptrs[bgcell] += 1
    end
    rewind_ptrs!(bgcell_to_cells_ptrs)
    Table(bgcell_to_cells_data,bgcell_to_cells_ptrs)
  end
  bgcell_to_bface = compress_to_bgcells(Γ,Ωᵇ)
  ccell_to_bface = lazy_map(Reindex(bgcell_to_bface),ccell_to_bgcell)

  # Combine cell and bface contributions
  c1 = array_cache(bface_mat)
  c2 = array_cache(bface_vec)
  cmat, cvec = unpair_arrays(ccell_to_matvec)
  cmat = lazy_map(cmat,ccell_to_bface) do mat1,bfaces
    mat = copy(mat1)
    for bface in bfaces
      mat2 = getindex!(c1,bface_mat,bface) # This is thread-unsafe (easy to fix though)
      mat .+= mat2
    end
    mat
  end
  cvec = lazy_map(cvec,ccell_to_bface) do vec1,bfaces
    vec = copy(vec1)
    for bface in bfaces
      vec2 = getindex!(c2,bface_vec,bface) # This is thread-unsafe (easy to fix though)
      vec .+= vec2
    end
    vec
  end
  cell_matvec = pair_arrays(cmat,cvec)
  cell_matvec = attach_constraints_cols(U,cell_matvec,ccell_to_bgcell)
  cell_matvec = attach_constraints_rows(V,cell_matvec,ccell_to_bgcell)

  # Find dofs at the interface
  # Here, dofs touched by more than one cell
  cell_dofs = lazy_map(Reindex(get_cell_dof_ids(V)),ccell_to_bgcell)
  dof_to_ncells = zeros(Int32,num_free_dofs(V))
  cache = array_cache(cell_dofs)
  for cell in 1:length(cell_dofs)
    dofs = getindex!(cache,cell_dofs,cell)
    for dof in dofs
      if dof > 0
        dof_to_ncells[dof] += 1
      end
    end
  end
  idof_to_dof = findall(ncells->ncells>1,dof_to_ncells)
  dof_to_idof = dof_to_ncells
  dof_to_idof .= -1
  dof_to_idof[idof_to_dof] = 1:length(idof_to_dof)

  # Local dofs at interface
  cell_ildof_ldof = lazy_map(cell_dofs) do dofs
    ildof_ldof = findall(dof->dof>0&&dof_to_idof[dof]>0,dofs)
    ildof_ldof
  end

  # Local dofs at cells
  cell_cldof_ldof = lazy_map(cell_dofs) do dofs
    cldof_ldof = findall(dof->dof>0&&dof_to_idof[dof]<0,dofs)
    cldof_ldof
  end

  # New local to global dof map
  cell_idofs = lazy_map(cell_dofs,cell_ildof_ldof) do dofs, ildof_ldof
    idofs = dof_to_idof[dofs[ildof_ldof]]
    idofs
  end

  # Static condensation
  cell_imatvec = lazy_map(
    cell_matvec,cell_ildof_ldof,cell_cldof_ldof) do matvec, ildof_ldof, cldof_ldof

    A, b = matvec
    Acc = A[cldof_ldof,cldof_ldof]
    Aci = A[cldof_ldof,ildof_ldof]
    Aic = A[ildof_ldof,cldof_ldof]
    Aii = A[ildof_ldof,ildof_ldof]
    bc = b[cldof_ldof]
    bi = b[ildof_ldof]
    Sii = Aii - Aic*(Acc\Aci)
    fi = bi - Aic*(Acc\bc)
    Sii, fi
  end

  # Assembly
  matvecdata =  ([cell_imatvec],[cell_idofs],[cell_idofs])
  matdata = ([],[],[])
  vecdata = ([],[])
  data = (matvecdata, matdata, vecdata)
  rows = first(axes(idof_to_dof))
  cols = rows
  assem = GenericSparseMatrixAssembler(
    SparseMatrixCSC{Float64,Int},
    Vector{Float64},
    rows,
    cols,
    DefaultAssemblyStrategy())
  S,f = assemble_matrix_and_vector(assem,data)

  # Solve
  x = S\f

  # Post-process
  y = zeros(num_free_dofs(V))
  y[idof_to_dof] = x
  arrays = (cell_matvec,cell_dofs,cell_ildof_ldof,cell_cldof_ldof)
  caches = map(array_cache,arrays)
  for cell in 1:length(cell_dofs)
    matvec,dofs,ildof_ldof,cldof_ldof = map((c,a)->getindex!(c,a,cell),caches,arrays)
    A, b = matvec
    Acc = A[cldof_ldof,cldof_ldof]
    Aci = A[cldof_ldof,ildof_ldof]
    bc = b[cldof_ldof]
    xi = x[dof_to_idof[dofs[ildof_ldof]]]
    xc = Acc\(bc - Aci*xi)
    y[dofs[cldof_ldof]] = xc
  end

  uh = FEFunction(V,y)
  writevtk(Ω,"trian",cellfields=["uh"=>uh])

  e = u - uh

  l2(u) = sqrt(sum( ∫( u*u )dΩ ))
  h1(u) = sqrt(sum( ∫( u*u + ∇(u)⋅∇(u) )dΩ ))

  el2 = l2(e)
  eh1 = h1(e)
  ul2 = l2(uh)
  uh1 = h1(uh)

  @test el2/ul2 < 1.e-8
  @test eh1/uh1 < 1.e-7

end # module
