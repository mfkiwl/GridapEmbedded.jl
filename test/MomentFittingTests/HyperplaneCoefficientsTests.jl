# module HyperplaneCoefficientsTests

using Gridap
using Gridap.Arrays
using GridapEmbedded.MomentFitting

facet_to_points_data = Int32[1,2,2,3,3,1,2,3,3,4,4,2]
facet_to_points_ptrs = Int32[1,3,5,7,9,11,13]
facet_to_points = Table(facet_to_points_data,facet_to_points_ptrs)

point_to_coords = VectorValue{2,Float64}[(0.0,0.0), (1.0,0.0), (0.0,1.0), (1.0,1.0)]

facet_to_normal = VectorValue{2,Float64}[(0.0,-1.0), (1/sqrt(2),1/sqrt(2)), (-1.0,0.0), (-1/sqrt(2),-1/sqrt(2)), (0.0,1.0), (1.0,0.0) ]

coeffs = compute_hyperplane_coefficients(facet_to_points,point_to_coords,facet_to_normal)

# coeffs = compute_hyperplane_coefficients(facet_to_points,point_to_coords)

# end #module
