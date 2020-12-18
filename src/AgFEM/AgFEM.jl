module AgFEM

using LightGraphs
using LinearAlgebra

using Gridap
using Gridap.Helpers
using Gridap.Arrays
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.CellData
using Gridap.FESpaces

using GridapEmbedded.CSG
using GridapEmbedded.Interfaces

export aggregate
export color_aggregates
export AggregateAllCutCells
export compute_aggregate_bboxes
export AgFEMSpace

include("CellAggregation.jl")

include("AggregateBoundingBoxes.jl")

include("AgFEMSpaces.jl")

end # module
