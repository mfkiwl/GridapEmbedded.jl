@inline function compute_hyperplane_coeffs(
                      facet_to_points::Table{Int32,Vector{Int32},Vector{Int32}},
                      point_to_coords::Vector{Point{D,T}},
                      facet_to_normal::AbstractArray) where {D,T}
  points = lazy_map(Reindex(facet_to_points.data),facet_to_points.ptrs[1:end-1])
  coords = lazy_map(Reindex(point_to_coords),points)
  lazy_map(⋅,facet_to_normal,coords)
end

@inline function compute_hyperplane_coeffs(
                      facet_to_points::Table{Int32,Vector{Int32},Vector{Int32}},
                      point_to_coords::Vector{Point{D,T}},
                      facet_to_normal::AbstractArray,
                      face_to_parent_face::AbstractArray) where {D,T}
  ptrs = lazy_map(Reindex(facet_to_points.ptrs),face_to_parent_face)
  points = lazy_map(Reindex(facet_to_points.data),ptrs)
  coords = lazy_map(Reindex(point_to_coords),points)
  lazy_map(⋅,facet_to_normal,coords)
end

function compute_hyperplane_coeffs(trian)
  @abstractmethod
end

function compute_hyperplane_coeffs(trian::SubFacetTriangulation)
  c = compute_hyperplane_coeffs(trian.subfacets)
  CellField(c,trian)
end

function compute_hyperplane_coeffs(sfd::SubFacetData)
  compute_hyperplane_coeffs(sfd.facet_to_points,
                            sfd.point_to_coords,
                            sfd.facet_to_normal)
end

function compute_hyperplane_coeffs(trian::SkeletonTriangulation)
  SkeletonPair(compute_hyperplane_coeffs(trian.⁺),
               compute_hyperplane_coeffs(trian.⁻))
end

function compute_hyperplane_coeffs(trian::SkeletonTriangulation,sym::Symbol)
  if sym in (:⁺,:plus)
    compute_hyperplane_coeffs(trian.⁺)
  elseif sym in (:⁻, :minus)
    compute_hyperplane_coeffs(trian.⁻)
  else
    @unreachable """\n
    It is not possible to use the given symbol on a SkeletonTriangulation.
    Make sure that you are specifying which of the two possible traces,
    either plus (aka ⁺) or minus (aka ⁻) you want to use.
    """
  end
end

function compute_hyperplane_coeffs(trian::SubFacetBoundaryTriangulation)
  c = compute_hyperplane_coeffs(trian.subfacets,trian.cell_normals)
  GenericCellField(c,trian,ReferenceDomain())
end

function compute_hyperplane_coeffs(scd::SubCellData,cell_normals::AbstractArray)
  compute_hyperplane_coeffs(scd.cell_to_points,
                            scd.point_to_coords,
                            cell_normals)
end

function compute_hyperplane_coeffs(trian::BoundaryTriangulation)
  facet_to_normal = get_facet_normal(trian)
  c = compute_hyperplane_coeffs(trian.face_trian,facet_to_normal)
  GenericCellField(c,trian,ReferenceDomain())
end

function compute_hyperplane_coeffs(trian::RestrictedTriangulation,
                                   cell_to_normal::AbstractArray)
  compute_hyperplane_coeffs(trian.parent_trian,
                            cell_to_normal,
                            trian.cell_to_parent_cell)
end

function compute_hyperplane_coeffs(trian::UnstructuredGrid,
                                   cell_to_normal::AbstractArray,
                                   cell_to_parent_cell::AbstractArray)
  compute_hyperplane_coeffs(trian.cell_node_ids,
                            trian.node_coordinates,
                            cell_to_normal,
                            cell_to_parent_cell)
end
