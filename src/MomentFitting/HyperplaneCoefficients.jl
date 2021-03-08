minor_indices(D::Int) = minor_indices(Val{D}())
minor_indices(::Val{D}) where {D} = @abstractmethod
minor_indices(::Val{2}) = [[2,3],[1,3],[1,2]]
minor_indices(::Val{3}) = [[1,2,3],[1,2,4],[1,3,4],[2,3,4]]

minor_sign(i::Int) = isodd(i) ? 1 : -1

function compute_hyperplane_coeffs(
  facet_to_points::Table{Int32,Vector{Int32},Vector{Int32}},point_to_coords::Vector{Point{D,T}}) where {D,T}
  coeffs = zeros(T,length(facet_to_points))
  c = array_cache(facet_to_points)
  fc = ones(T,length(c),D+1)
  m = minor_indices(D)
  for i = 1:length(facet_to_points)
    getindex!(c,facet_to_points,i)
    for j in 1:length(c)
      fc[j,1:D] .= point_to_coords[c][j].data
    end
    a = 0.0
    for j in 1:D
      a += det(view(fc,:,m[j]))^2
    end
    b = -minor_sign(D+1)*det(view(fc,:,m[D+1]))
    # Use precomputed exterior normal to orient hyperplane
    coeffs[i] = b / sqrt(a)
  end
  coeffs
end

@inline function compute_hyperplane_coeffs(
  facet_to_points::Table{Int32,Vector{Int32},Vector{Int32}},point_to_coords::Vector{Point{D,T}},
  facet_to_normal::AbstractArray) where {D,T}
  points = lazy_map(Reindex(facet_to_points.data),facet_to_points.ptrs[1:end-1])
  coords = lazy_map(Reindex(point_to_coords),points)
  lazy_map(⋅,facet_to_normal,coords)
end

@inline function compute_hyperplane_coeffs(
  facet_to_points::Table{Int32,Vector{Int32},Vector{Int32}},point_to_coords::Vector{Point{D,T}},
  facet_to_normal::AbstractArray,
  face_to_parent_face::AbstractArray) where {D,T}
  ptrs = lazy_map(Reindex(facet_to_points.ptrs),face_to_parent_face)
  points = lazy_map(Reindex(facet_to_points.data),ptrs)
  coords = lazy_map(Reindex(point_to_coords),points)
  lazy_map(⋅,facet_to_normal,coords)
end

function compute_hyperplane_coeffs(trian::SubFacetTriangulation)
  compute_hyperplane_coeffs(trian.subfacets)
end

function compute_hyperplane_coeffs(sfd::SubFacetData)
  compute_hyperplane_coeffs(sfd.facet_to_points,
                            sfd.point_to_coords,
                            sfd.facet_to_normal)
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
  compute_hyperplane_coeffs(trian.subfacets,trian.cell_normals)
end

function compute_hyperplane_coeffs(scd::SubCellData,cell_normals::AbstractArray)
  compute_hyperplane_coeffs(scd.cell_to_points,
                            scd.point_to_coords,
                            cell_normals)
end

function compute_hyperplane_coeffs(trian::BoundaryTriangulation)
  facet_to_normal = get_facet_normal(trian)
  compute_hyperplane_coeffs(trian.face_trian,facet_to_normal)
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
