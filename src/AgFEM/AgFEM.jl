module AgFEM

using LightGraphs
using LinearAlgebra

using Gridap
using Gridap.Helpers
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Geometry: get_cell_to_parent_cell
using Gridap.ReferenceFEs
using Gridap.CellData
using Gridap.FESpaces
using Gridap.FESpaces: VoidBasis

using GridapEmbedded.CSG
using GridapEmbedded.Interfaces

export aggregate
export color_aggregates
export AggregateAllCutCells
export AgFEMSpace

include("CellAggregation.jl")

include("AgFEMSpaces.jl")

end # module
