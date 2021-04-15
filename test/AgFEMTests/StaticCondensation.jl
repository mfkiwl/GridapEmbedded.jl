  using Gridap
  using Gridap.Geometry
  using Gridap.CellData
  using Gridap.Arrays
  using Gridap.Algebra
  using SparseMatricesCSR
  using Gridap.FESpaces

  import Gridap.Arrays: return_cache, evaluate!

  function compute_contributions(U::FESpace,V::FESpace,
                                 a::Function,l::Function,
                                 Ω::Triangulation,Γ::Triangulation)
    du = get_cell_shapefuns_trial(U)
    dv = get_cell_shapefuns(V)
    cell_mat = a(du,dv)[Ω]
    cell_vec = l(dv)[Ω]
    bface_mat = a(du,dv)[Γ]
    bface_vec = l(dv)[Γ]
    pair_arrays(cell_mat,cell_vec), pair_arrays(bface_mat,bface_vec)
  end

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

  struct CellArrayMap{A} <: Map
    array::A
  end

  function return_cache(a::CellArrayMap,array,bfaces)
    ca = copy(array)
    c = array_cache(a.array)
    ca, c
  end

  function evaluate!(cache,a::CellArrayMap,array,bfaces)
    ca, c = cache
    copyto!(ca,array)
    for bface in bfaces
      ab = getindex!(c,a.array,bface)
      ca .+= ab
    end
    ca
  end

  function combine_cell_and_bface_contribs(bgtrian::Triangulation,
                                           trian::Triangulation,
                                           ebtrian::Triangulation,
                                           cell_matvec::AbstractArray,
                                           bface_matvec::AbstractArray)
    cell_to_bgcell = get_cell_to_bgcell(trian)
    ccell_to_first_cell = compress(cell_to_bgcell)
    ccell_to_matvec, ccell_to_bgcell = compress(cell_matvec,
                                                cell_to_bgcell,
                                                ccell_to_first_cell)
    bgcell_to_bface = compress_to_bgcells(ebtrian,bgtrian)
    ccell_to_bface = lazy_map(Reindex(bgcell_to_bface),ccell_to_bgcell)
    cmat, cvec = unpair_arrays(ccell_to_matvec)
    bmat, bvec = unpair_arrays(bface_matvec)
    cmat = lazy_map(CellArrayMap(bmat),cmat,ccell_to_bface)
    cvec = lazy_map(CellArrayMap(bvec),cvec,ccell_to_bface)
    pair_arrays(cmat,cvec), ccell_to_bgcell
  end

  function compute_itfc_dofs(cell_dofs::AbstractArray,V::FESpace)
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
    idof_to_dof, dof_to_idof
  end

  function compute_itfc_cell_dof_arrays(cell_dofs::AbstractArray,
                                        dof_to_idof::AbstractArray)
    cell_ildof_ldof = lazy_map(cell_dofs) do dofs
      ildof_ldof = findall(dof->dof>0&&dof_to_idof[dof]>0,dofs)
      ildof_ldof
    end
    cell_cldof_ldof = lazy_map(cell_dofs) do dofs
      cldof_ldof = findall(dof->dof>0&&dof_to_idof[dof]<0,dofs)
      cldof_ldof
    end
    cell_idofs = lazy_map(cell_dofs,cell_ildof_ldof) do dofs, ildof_ldof
      idofs = dof_to_idof[dofs[ildof_ldof]]
      idofs
    end
    cell_ildof_ldof, cell_cldof_ldof, cell_idofs
  end

  struct StaticCondensationMap <: Map
    dict::Dict{Any,Any}
  end

  function _return_cache(dict::Dict,key::Tuple{Int64,Int64})
    id = objectid(key)
    if ! haskey(dict,id)
      dict[id] = zeros(Float64,key[1],key[2])
    end
    dict[id]
  end

  function return_cache(a::StaticCondensationMap,matvec,ildof_ldof,cldof_ldof)
    cacc = _return_cache(a.dict,(length(cldof_ldof),length(cldof_ldof)))
    caci = _return_cache(a.dict,(length(cldof_ldof),length(ildof_ldof)))
    caic = _return_cache(a.dict,(length(ildof_ldof),length(cldof_ldof)))
    caii = _return_cache(a.dict,(length(ildof_ldof),length(ildof_ldof)))
    cacc, caci, caic, caii
  end

  function evaluate!(cache,a::StaticCondensationMap,matvec,ildof_ldof,cldof_ldof)
    cacc, caci, caic, caii = return_cache(a::StaticCondensationMap,matvec,ildof_ldof,cldof_ldof)
    A, b = matvec
    if length(cldof_ldof) > 0
      copyto!(cacc,view(A,cldof_ldof,cldof_ldof))
      copyto!(caci,view(A,cldof_ldof,ildof_ldof))
      copyto!(caic,view(A,ildof_ldof,cldof_ldof))
      copyto!(caii,view(A,ildof_ldof,ildof_ldof))
      cbc = b[cldof_ldof]
      cbi = b[ildof_ldof]
      Sii = caii - caic*(cacc\caci)
      fi = cbi - caic*(cacc\cbc)
      Sii, fi
    else
      A, b
    end
  end


  function assemble_schur_system(cell_imatvec::AbstractArray,
                                 cell_idofs::AbstractArray,
                                 idof_to_dof::AbstractArray)
    matvecdata =  ([cell_imatvec],[cell_idofs],[cell_idofs])
    matdata = ([],[],[])
    vecdata = ([],[])
    data = (matvecdata, matdata, vecdata)
    rows = first(axes(idof_to_dof))
    cols = rows
    # assem = GenericSparseMatrixAssembler(
    #   SparseMatrixCSC{Float64,Int},
    #   Vector{Float64},
    #   rows,
    #   cols,
    #   DefaultAssemblyStrategy())
    assem = GenericSparseMatrixAssembler(
      SparseMatrixCSR{0,PetscReal,PetscInt},
      Vector{Float64},
      rows,
      cols,
      DefaultAssemblyStrategy())
    S,f = assemble_matrix_and_vector(assem,data)
  end

  function apply_interior_correction(x::Array{Float64,1},
                                     V::FESpace,
                                     idof_to_dof::AbstractArray,
                                     dof_to_idof::AbstractArray,
                                     cell_matvec::AbstractArray,
                                     cell_dofs::AbstractArray,
                                     cell_ildof_ldof::AbstractArray,
                                     cell_cldof_ldof::AbstractArray,
                                     dict::Dict)
    y = zeros(num_free_dofs(V))
    y[idof_to_dof] = x
    arrays = (cell_matvec,cell_dofs,cell_ildof_ldof,cell_cldof_ldof)
    caches = map(array_cache,arrays)
    for cell in 1:length(cell_dofs)
      matvec,dofs,ildof_ldof,cldof_ldof = map((c,a)->getindex!(c,a,cell),caches,arrays)
      if length(cldof_ldof) > 0
        A, b = matvec
        Acc = _return_cache(dict,(length(cldof_ldof),length(cldof_ldof)))
        Aci = _return_cache(dict,(length(cldof_ldof),length(ildof_ldof)))
        copyto!(Acc,view(A,cldof_ldof,cldof_ldof))
        copyto!(Aci,view(A,cldof_ldof,ildof_ldof))
        bc = b[cldof_ldof]
        xi = x[dof_to_idof[dofs[ildof_ldof]]]
        xc = Acc\(bc - Aci*xi)
        y[dofs[cldof_ldof]] = xc
      end
    end
    y
  end
