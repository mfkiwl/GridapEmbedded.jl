module MomentFitting

using Gridap
using Gridap.Helpers
using Gridap.Arrays
using Gridap.Fields: linear_combination

using Gridap.CellData

using Gridap.Geometry
using Gridap.Geometry: get_cell_to_parent_cell
using Gridap.Geometry: FaceToCellGlue

using Gridap.Polynomials

using Gridap.ReferenceFEs: compute_nodes
using Gridap.ReferenceFEs: LagrangianDofBasis

using GridapEmbedded.Interfaces
using GridapEmbedded.Interfaces: SubFacetTriangulation
using GridapEmbedded.Interfaces: SubFacetBoundaryTriangulation

using GridapEmbedded.LevelSetCutters: _check_and_get_polytope

import LinearAlgebra: I
using FillArrays

export compute_cell_moments

include("HyperplaneCoefficients.jl")

include("CutCellMoments.jl")

end #module
