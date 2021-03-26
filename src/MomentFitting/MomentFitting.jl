module MomentFitting

using Gridap
using Gridap.Helpers
using Gridap.Arrays
using Gridap.Fields: linear_combination, IntegrationMap

using Gridap.CellData

using Gridap.Geometry
using Gridap.Geometry: get_reffes
using Gridap.Geometry: get_cell_to_parent_cell
using Gridap.Geometry: FaceToCellGlue

using Gridap.Polynomials

using Gridap.ReferenceFEs: get_polytope
using Gridap.ReferenceFEs: compute_nodes
using Gridap.ReferenceFEs: LagrangianDofBasis

using Gridap.Arrays: CompressedArray

using Gridap.Integration: Quadrature
using Gridap.Integration: GenericQuadrature

using GridapEmbedded.Interfaces
using GridapEmbedded.Interfaces: SubFacetTriangulation
using GridapEmbedded.Interfaces: SubFacetBoundaryTriangulation

import GridapEmbedded.Interfaces: compute_bgcell_to_inoutcut

using GridapEmbedded.LevelSetCutters: _check_and_get_polytope

import LinearAlgebra: I, pinv, cond
using FillArrays

export MomentFittingMeasures
export MomentFittingQuad

include("HyperplaneCoefficients.jl")

include("CutCellMoments.jl")

end #module
