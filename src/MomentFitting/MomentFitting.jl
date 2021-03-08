module MomentFitting

using Gridap
using Gridap.Helpers
using Gridap.Arrays

using Gridap.CellData: DomainContribution

using GridapEmbedded.Interfaces: SubFacetTriangulation, SubFacetData
using GridapEmbedded.Interfaces: SubFacetBoundaryTriangulation, SubCellData

using Gridap.Geometry: UnstructuredGrid, get_facet_normal
using Gridap.Geometry: RestrictedDiscreteModel, get_cell_to_parent_cell

using FillArrays

export compute_hyperplane_coeffs
export compute_cell_moments

include("HyperplaneCoefficients.jl")

include("CutCellMoments.jl")

end #module
