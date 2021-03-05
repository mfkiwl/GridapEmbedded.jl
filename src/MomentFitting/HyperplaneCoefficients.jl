minor_indices(D::Int) = minor_indices(Val{D}())
minor_indices(::Val{D}) where {D} = @abstractmethod
minor_indices(::Val{2}) = [[2,3],[1,3],[1,2]]
minor_indices(::Val{3}) = [[1,2,3],[1,2,4],[1,3,4],[2,3,4]]

minor_sign(i::Int) = isodd(i) ? 1 : -1

function compute_hyperplane_coefficients(
  facet_to_points::Table{Int32,Vector{Int32},Vector{Int32}},point_to_coords::Vector{Point{D,T}},
  facet_to_normal::Vector{VectorValue{D,T}}) where {D,T}
  points = lazy_map(Reindex(facet_to_points.data),facet_to_points.ptrs[1:end-1])
  coords = lazy_map(Reindex(point_to_coords),points)
  lazy_map(⋅,facet_to_normal,coords)
end

function compute_hyperplane_coefficients(
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
