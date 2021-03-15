struct CutCellMoments
  data::Vector{Vector{Float64}}
  bgcell_to_cut_cell::Vector{Int32}
end

function CutCellMoments(model::RestrictedDiscreteModel,
                        facet_moments::DomainContribution)
  fi = [ testitem(array) for (trian,array) in facet_moments.dict ]
  li = map(length,fi)
  @assert all(li .== first(li))
  cell_to_parent_cell = get_cell_to_parent_cell(model)
  data = [ zero(first(fi)) for i in 1:length(cell_to_parent_cell) ]
  bgcell_to_cut_cell = zeros(Int32,num_cells(get_parent_model(model)))
  bgcell_to_cut_cell[cell_to_parent_cell] .= 1:length(cell_to_parent_cell)
  CutCellMoments(data,bgcell_to_cut_cell)
end

function MomentFittingQuad(active_mesh::RestrictedTriangulation,
                           cut::EmbeddedDiscretization,
                           degree::Int)

  acell_to_point_vals, acell_to_weight_vals = compute_cell_moments(cut,degree)
  acell_to_weight_vals = collect(get_array(acell_to_weight_vals))

  bgcell_to_inoutcut = compute_bgcell_to_inoutcut(cut)
  acell_to_bgcell = active_mesh.cell_to_parent_cell
  acell_to_inoutcut = lazy_map(Reindex(bgcell_to_inoutcut),acell_to_bgcell)
  acell_to_point_ptrs = lazy_map(i->(i == CUT ? 1 : 2),acell_to_inoutcut)

  quad = map(r->Quadrature(get_polytope(r),2*degree),get_reffes(active_mesh))
  @assert length(quad) == 1
  acell_to_point_vals = [acell_to_point_vals,get_coordinates(quad[1])]

  push!(acell_to_weight_vals,get_weights(quad[1]))

  acell_to_is_cut = findall(lazy_map(i->(i == CUT),acell_to_inoutcut))
  num_quads = length(acell_to_weight_vals)
  acell_to_weight_ptrs = map(i->(i == IN ? num_quads : 0),acell_to_inoutcut)
  acell_to_weight_ptrs[acell_to_is_cut] .= 1:length(acell_to_is_cut)

  acell_to_point = CompressedArray(acell_to_point_vals,acell_to_point_ptrs)
  acell_to_weight = CompressedArray(acell_to_weight_vals,acell_to_weight_ptrs)

  acell_to_quad = [ GenericQuadrature(acell_to_point[i],acell_to_weight[i]) for i in 1:num_quads ]
  CellQuadrature(acell_to_quad,acell_to_point,acell_to_weight,active_mesh,ReferenceDomain())

end

function compute_cell_moments(cut::EmbeddedDiscretization{D,T},
                              degree::Int) where{D,T}
  bgtrian = Triangulation(cut.bgmodel)
  b = MonomialBasis{D}(T,degree)
  v = get_monomial_cell_field(b,bgtrian)
  cut_bgmodel = DiscreteModel(cut,cut.geo,CUT)
  mon_contribs = compute_monomial_domain_contribution(cut,degree,v)
  mon_moments = compute_monomial_cut_cell_moments(cut_bgmodel,mon_contribs,b)
  lag_nodes, lag_to_mon = get_nodes_and_change_of_basis(cut_bgmodel,cut,b,degree)
  lag_moments = lazy_map(*,lag_to_mon,mon_moments)
  lag_moments = map_to_ref_space!(lag_moments,lag_nodes,cut_bgmodel)
  lag_nodes, lag_moments
end

function compute_monomial_domain_contribution(cut::EmbeddedDiscretization{D,T},
                                              degree::Int,
                                              v::CellField) where {D,T}

  Γᵉ = EmbeddedBoundary(cut)
  Λ  = GhostSkeleton(cut)
  cutf = cut_facets(cut.bgmodel,cut.geo)
  Γᶠ = SkeletonTriangulation(cutf,Λ,cut.geo,CUTIN)
  Γᵇ = SkeletonTriangulation(cutf,Λ,cut.geo,IN)
  Γᵒ = BoundaryTriangulation(cutf,CUTIN)
  Λ  = BoundaryTriangulation(cut.bgmodel)
  Γᵖ = BoundaryTriangulation(cutf,Λ,cut.geo,IN)

  cutdeg = 2*D*degree
  dΓᵉ = Measure(Γᵉ,cutdeg)
  dΓᶠ = SkeletonPair(Measure(Γᶠ.⁺,cutdeg),Measure(Γᶠ.⁻,cutdeg))
  dΓᵇ = SkeletonPair(Measure(Γᵇ.⁺,cutdeg),Measure(Γᵇ.⁻,cutdeg))
  dΓᵒ = Measure(Γᵒ,cutdeg)
  dΓᵖ = Measure(Γᵖ,cutdeg)

  cᵉ = compute_hyperplane_coeffs(Γᵉ)
  cᶠ = compute_hyperplane_coeffs(Γᶠ)
  cᵇ = compute_hyperplane_coeffs(Γᵇ)
  cᵒ = compute_hyperplane_coeffs(Γᵒ)
  cᵖ = compute_hyperplane_coeffs(Γᵖ)

  @check num_cells(Γᵉ) > 0
  J = ∫(cᵉ*v)*dΓᵉ +
      ∫(cᶠ.⁺*v)*dΓᶠ.⁺ + ∫(cᶠ.⁻*v)*dΓᶠ.⁻
  if num_cells(Γᵇ) > 0
    J += ∫(cᵇ.⁺*v)*dΓᵇ.⁺ + ∫(cᵇ.⁻*v)*dΓᵇ.⁻
  end
  if num_cells(Γᵒ) > 0
    J += ∫(cᵒ*v)*dΓᵒ
  end
  if num_cells(Γᵖ) > 0
    J += ∫(cᵖ*v)*dΓᵖ
  end
  J

end

function compute_monomial_cut_cell_moments(model::RestrictedDiscreteModel,
                                           facet_moments::DomainContribution,
                                           b::MonomialBasis{D,T}) where {D,T}
  cut_cell_to_moments = CutCellMoments(model,facet_moments)
  for (trian,array) in facet_moments.dict
    add_facet_moments!(cut_cell_to_moments,trian,array)
  end
  o = get_terms_degrees(b)
  q = 1 ./ ( D .+ o )
  [ q .* d for d in cut_cell_to_moments.data ]
end

function add_facet_moments!(ccm::CutCellMoments,trian,array::AbstractArray)
  @abstractmethod
end

function add_facet_moments!(ccm::CutCellMoments,
                            trian::SubFacetTriangulation,
                            array::AbstractArray)
  add_facet_moments!(ccm,trian.subfacets,array)
end

function add_facet_moments!(ccm::CutCellMoments,
                            sfd::SubFacetData,
                            array::AbstractArray)
  facet_to_cut_cell = lazy_map(Reindex(ccm.bgcell_to_cut_cell),sfd.facet_to_bgcell)
  for i = 1:length(facet_to_cut_cell)
    ccm.data[facet_to_cut_cell[i]] += array[i]
  end
end

function add_facet_moments!(ccm::CutCellMoments,
                            trian::SubFacetBoundaryTriangulation,
                            array::AbstractArray)
  if length(trian.subfacet_to_facet) > 0
    subfacet_to_bgcell = lazy_map(Reindex(trian.facets.glue.face_to_cell),trian.subfacet_to_facet)
    subfacet_to_cut_cell = lazy_map(Reindex(ccm.bgcell_to_cut_cell),subfacet_to_bgcell)
    l = length(subfacet_to_cut_cell)
    for i = 1:l
      ccm.data[subfacet_to_cut_cell[i]] += array[i]
    end
  else
    add_facet_moments!(ccm,trian.facets,array)
  end
end

function add_facet_moments!(ccm::CutCellMoments,
                            trian::BoundaryTriangulation,
                            array::AbstractArray)
  add_facet_moments!(ccm,trian.glue,array)
end

function add_facet_moments!(ccm::CutCellMoments,
                            glue::FaceToCellGlue,
                            array::AbstractArray)
  facet_to_cut_cell = lazy_map(Reindex(ccm.bgcell_to_cut_cell),glue.face_to_cell)
  cell_to_is_cut = findall(lazy_map(i->(i>0),facet_to_cut_cell))
  facet_to_cut_cell = lazy_map(Reindex(facet_to_cut_cell),cell_to_is_cut)
  l = length(facet_to_cut_cell)
  for i = 1:l
    ccm.data[facet_to_cut_cell[i]] += array[cell_to_is_cut[i]]
  end
end

function get_nodes_and_change_of_basis(model::RestrictedDiscreteModel,
                                       cut::EmbeddedDiscretization{D,T},
                                       b::MonomialBasis{D,T},
                                       degree::Int) where {D,T}
  p = check_and_get_polytope(cut)
  orders = tfill(degree,Val{D}())
  nodes, _ = compute_nodes(p,orders)
  coords = get_cell_coordinates(model)
  dofs = [ LagrangianDofBasis(T,cs) for cs in coords ]
  change = [ transpose((inv(evaluate(ds,b)))) for ds in dofs ]
  nodes, change
end

function compute_bgcell_to_inoutcut(cut::EmbeddedDiscretization)
  compute_bgcell_to_inoutcut(cut,cut.geo)
end

function map_to_ref_space!(moments::AbstractArray,
                           nodes::Vector{<:Point},
                           model::RestrictedDiscreteModel)
  cell_map = get_cell_map(model)
  cell_Jt = lazy_map(∇,cell_map)
  cell_detJt = lazy_map(Operation(det),cell_Jt)
  cell_nodes = Fill(nodes,num_cells(model))
  detJt = lazy_map(evaluate,cell_detJt,cell_nodes)
  moments = lazy_map(Broadcasting(/),moments,detJt)
end

@inline function check_and_get_polytope(cut::EmbeddedDiscretization)
  _check_and_get_polytope(cut.bgmodel.grid)
end

function get_monomial_cell_field(b::MonomialBasis{D,T},
                                 trian::Triangulation) where {D,T}
  i = Matrix{eltype(T)}(I,length(b),length(b))
  l = linear_combination(i,b)
  m = Fill(l,num_cells(trian))
  GenericCellField(m,trian,PhysicalDomain())
end

@inline function get_terms_degrees(b::MonomialBasis)
  [ _get_terms_degrees(c) for c in b.terms ]
end

function _get_terms_degrees(c::CartesianIndex)
  d = 0
  for i in 1:length(c)
    d += (c[i]-1)
  end
  d
end
