
function AgFEMSpace(f::SingleFieldFESpace,cell_to_cellin::AbstractVector,g::SingleFieldFESpace=f)
  AgFEMSpace(f,cell_to_cellin,get_cell_shapefuns(g),get_cell_dof_basis(g))
end

# Note: cell is in fact bgcell in this function since f will usually be an ExtendedFESpace
function AgFEMSpace(
  f::SingleFieldFESpace,
  cell_to_cellin::AbstractVector,
  cell_shapefuns_g::CellField,
  cell_dof_basis_g::CellDof)

  # Prepare maps between different cell ids
  cell_to_isactive = lazy_map(i->(i>0),cell_to_cellin)
  acell_to_cell = findall( cell_to_isactive  )
  acell_to_cellin = cell_to_cellin[acell_to_cell]
  cell_to_acell = zeros(Int32,length(cell_to_cellin))
  cell_to_acell[cell_to_isactive] .= 1:length(acell_to_cell)

  # Triangulation made of active cells
  trian = get_triangulation(f)
  trian_a = Triangulation(trian,acell_to_cell)

  # Celldata Defined on trian
  dofs_f = get_cell_dof_basis(f)
  shfns_f = get_cell_shapefuns(f)
  dofs_g = cell_dof_basis_g
  shfns_g = cell_shapefuns_g

  # Celldata Defined on trian_a
  dofs_f_a = change_domain(dofs_f,trian_a,DomainStyle(dofs_f))
  shfns_g_a = lazy_map(Reindex(get_data(shfns_g)),acell_to_cellin)
  shfns_g_a = GenericCellField(shfns_g_a,trian_a,DomainStyle(dofs_f_a))

  # LAGRANGIAN
  # dofs_f_a_data = get_data(dofs_f_a)
  # dofs_f_a_bases = map(i->i.basis,dofs_f_a_data)
  # dofs_f_a_isvoid = map(i->i.isvoid,dofs_f_a_data)

  # cell_ref_nodes = map(i->i.nodes,dofs_f_a_bases)
  # cell_map = get_cell_map(trian_a)
  # cell_phys_nodes = lazy_map(evaluate,cell_map,cell_ref_nodes)

  # cell_map = get_cell_map(trian)
  # acell_map = lazy_map(Reindex(cell_map),acell_to_cellin)
  # acell_invmap = lazy_map(inverse_map,acell_map)
  # acell_ref_nodes = lazy_map(evaluate,acell_invmap,cell_phys_nodes)

  # adofs_f_a_bases = lazy_map(LagrangianDofBasis,dofs_f_a_bases,acell_ref_nodes)
  # adofs_f_a_data = lazy_map(VoidBasis,adofs_f_a_bases,dofs_f_a_isvoid)
  # adofs_f_a = CellDof(adofs_f_a_data,trian_a,DomainStyle(dofs_f_a))
  # END LAGRANGIAN

  # MODAL
  dofs_f_a_data = get_data(dofs_f_a)
  dofs_f_a_bases = map(i->i.basis,dofs_f_a_data)
  lag_dofs_f_a_bases = map(i->i.basis.dof_basis,dofs_f_a_data)
  dofs_f_a_isvoid = map(i->i.isvoid,dofs_f_a_data)

  cell_ref_nodes = map(i->i.nodes,lag_dofs_f_a_bases)
  cell_map = get_cell_map(trian_a)
  cell_phys_nodes = lazy_map(evaluate,cell_map,cell_ref_nodes)

  cell_map = get_cell_map(trian)
  acell_map = lazy_map(Reindex(cell_map),acell_to_cellin)
  acell_invmap = lazy_map(inverse_map,acell_map)
  acell_ref_nodes = lazy_map(evaluate,acell_invmap,cell_phys_nodes)

  lag_adofs_f_a_bases = lazy_map(LagrangianDofBasis,lag_dofs_f_a_bases,acell_ref_nodes)
  adofs_f_a_bases = lazy_map(linear_combination,dofs_f_a_bases,lag_adofs_f_a_bases)
  adofs_f_a_data = lazy_map(VoidBasis,adofs_f_a_bases,dofs_f_a_isvoid)
  adofs_f_a = CellDof(adofs_f_a_data,trian_a,DomainStyle(dofs_f_a))
  # END MODAL

  # dofs_f_a = change_domain(dofs_f,trian_a,DomainStyle(dofs_f))
  # cell_phys_shapefuns_g = get_array(change_domain(shfns_g,PhysicalDomain()))
  # acell_phys_shapefuns_g = lazy_map(Reindex(cell_phys_shapefuns_g),acell_to_cellin)
  # shfns_g_a = GenericCellField(acell_phys_shapefuns_g,trian_a,PhysicalDomain())

  # acell_to_coeffs = dofs_f_a(shfns_g_a)
  acell_to_coeffs = adofs_f_a(shfns_g_a)
  cell_to_proj = dofs_g(shfns_f)
  acell_to_proj = lazy_map(Reindex(cell_to_proj),acell_to_cellin)
  acell_to_dof_ids = lazy_map(Reindex(get_cell_dof_ids(f)),acell_to_cell)

  #acell_to_dofs = reindex(get_cell_dofs(f),acell_to_cell)
  #n_fdofs = 
  #acell_to_fbasis = reindex(get_cell_basis(f),acell_to_cellin)
  #acell_to_gbasis = reindex(cell_shapefuns_g,acell_to_cellin)
  #acell_to_dof_fbasis = reindex(get_cell_dof_basis(f),acell_to_cell)
  #acell_to_dof_gbasis = reindex(cell_dof_basis_g,acell_to_cellin)
  #@notimplementedif is_in_ref_space(acell_to_dof_fbasis)
  #@notimplementedif is_in_ref_space(acell_to_gbasis)
  #@assert RefStyle(acell_to_fbasis) == RefStyle(acell_to_dof_gbasis)
  #acell_to_coeffs = evaluate(acell_to_dof_fbasis,acell_to_gbasis)
  #acell_to_proj = evaluate(acell_to_dof_gbasis,acell_to_fbasis)

  aggdof_to_fdof, aggdof_to_dofs, aggdof_to_coeffs = _setup_agfem_constraints(
    num_free_dofs(f),
    acell_to_cellin,
    acell_to_cell,
    cell_to_acell,
    acell_to_dof_ids,
    acell_to_coeffs,
    acell_to_proj)

  FESpaceWithLinearConstraints(aggdof_to_fdof,aggdof_to_dofs,aggdof_to_coeffs,f)
end

function _setup_agfem_constraints(
  n_fdofs,
  acell_to_cellin,
  acell_to_cell,
  cell_to_acell,
  acell_to_dof_ids,
  acell_to_coeffs,
  acell_to_proj)

  n_acells = length(acell_to_cell)
  fdof_to_isagg = fill(true,n_fdofs)
  fdof_to_acell = zeros(Int32,n_fdofs)
  fdof_to_ldof = zeros(Int8,n_fdofs)
  cache = array_cache(acell_to_dof_ids)
  for acell in 1:n_acells
    iscut = acell_to_cell[acell] != acell_to_cellin[acell]
    dofs = getindex!(cache,acell_to_dof_ids,acell)
    cellin = acell_to_cellin[acell]
    for (ldof,dof) in enumerate(dofs)
      if dof > 0
        fdof = dof
        fdof_to_isagg[fdof] = iscut && fdof_to_isagg[fdof]
        fdof_to_acell[fdof] = acell
        fdof_to_ldof[fdof] = ldof
      end
    end
  end

  aggdof_to_fdof = findall(fdof_to_isagg)

  n_aggdofs = length(aggdof_to_fdof)
  aggdof_to_dofs_ptrs = zeros(Int32,n_aggdofs+1)

  for aggdof in 1:n_aggdofs
    fdof = aggdof_to_fdof[aggdof]
    acell = fdof_to_acell[fdof]
    cellin = acell_to_cellin[acell]
    acell = cell_to_acell[cellin]
    dofs = getindex!(cache,acell_to_dof_ids,acell)
    aggdof_to_dofs_ptrs[aggdof+1] = length(dofs)
  end

  length_to_ptrs!(aggdof_to_dofs_ptrs)
  ndata = aggdof_to_dofs_ptrs[end]-1
  aggdof_to_dofs_data = zeros(Int,ndata)

  for aggdof in 1:n_aggdofs
    fdof = aggdof_to_fdof[aggdof]
    acell = fdof_to_acell[fdof]
    cellin = acell_to_cellin[acell]
    acell = cell_to_acell[cellin]
    dofs = getindex!(cache,acell_to_dof_ids,acell)
    p = aggdof_to_dofs_ptrs[aggdof]-1
    for (i,dof) in enumerate(dofs)
      aggdof_to_dofs_data[p+i] = dof
    end
  end

  aggdof_to_dofs = Table(aggdof_to_dofs_data,aggdof_to_dofs_ptrs)

  cache2 = array_cache(acell_to_coeffs)
  cache3 = array_cache(acell_to_proj)

  T = eltype(eltype(acell_to_coeffs))
  z = zero(T)

  aggdof_to_coefs_data = zeros(T,ndata)
  for aggdof in 1:n_aggdofs
    fdof = aggdof_to_fdof[aggdof]
    acell = fdof_to_acell[fdof]
    coeffs = getindex!(cache2,acell_to_coeffs,acell)
    proj = getindex!(cache3,acell_to_proj,acell)
    ldof = fdof_to_ldof[fdof]
    p = aggdof_to_dofs_ptrs[aggdof]-1
    for b in 1:size(proj,2)
      coeff = z
      for c in 1:size(coeffs,2)
        coeff += coeffs[ldof,c]*proj[c,b]
      end
      aggdof_to_coefs_data[p+b] = coeff
    end
  end

  aggdof_to_coeffs = Table(aggdof_to_coefs_data,aggdof_to_dofs_ptrs)

  aggdof_to_fdof, aggdof_to_dofs, aggdof_to_coeffs
end

