module AggregateBoundingBoxesTests

using Gridap
using GridapEmbedded
using GridapEmbedded.AgFEM

const R = 0.4

geom = disk(R,x0=Point(0.5,0.5))
n = 5
partition = (n,n)
domain = (0,1,0,1)

# geom = sphere(R,x0=Point(0.5,0.5,0.5))
# n = 5
# partition = (n,n,n)
# domain = (0,1,0,1,0,1)

bgmodel = CartesianDiscreteModel(domain,partition)

cutdisc = cut(bgmodel,geom)

strategy = AggregateAllCutCells()

aggregates = aggregate(strategy,cutdisc,geom)

colors = color_aggregates(aggregates,bgmodel)

trian = Triangulation(bgmodel)

coords = collect(get_cell_coordinates(trian))

bboxes = compute_aggregate_bboxes(bgmodel,aggregates)

anchors = compute_anchor_nfaces(bgmodel)

dbboxes = compute_bbox_dfaces(bgmodel,bboxes)

sktrian = BoundaryTriangulation(bgmodel,tags=collect(1:9))

bbminsx = [ dbboxes[1][i][1].data[1] for i in 1:length(dbboxes[1]) ]
bbminsy = [ dbboxes[1][i][1].data[2] for i in 1:length(dbboxes[1]) ]

bbmaxsx = [ dbboxes[1][i][2].data[1] for i in 1:length(dbboxes[1]) ]
bbmaxsy = [ dbboxes[1][i][2].data[2] for i in 1:length(dbboxes[1]) ]

writevtk(sktrian,"sktrian",
         celldata=["aminsx"=>bbminsx,"aminsy"=>bbminsy,
                   "bmaxsx"=>bbmaxsx,"bmaxsy"=>bbmaxsy])

# sktrian = BoundaryTriangulation(bgmodel,tags=collect(1:27))

# bbminsx = [ dbboxes[1][i][1].data[1] for i in 1:length(dbboxes[2]) ]
# bbminsy = [ dbboxes[1][i][1].data[2] for i in 1:length(dbboxes[2]) ]
# bbminsz = [ dbboxes[1][i][1].data[3] for i in 1:length(dbboxes[2]) ]

# bbmaxsx = [ dbboxes[1][i][2].data[1] for i in 1:length(dbboxes[2]) ]
# bbmaxsy = [ dbboxes[1][i][2].data[2] for i in 1:length(dbboxes[2]) ]
# bbmaxsz = [ dbboxes[1][i][2].data[3] for i in 1:length(dbboxes[2]) ]

# writevtk(sktrian,"sktrian",
#          celldata=["aminsx"=>bbminsx,"aminsy"=>bbminsy,"aminsz"=>bbminsz,
#                    "bmaxsx"=>bbmaxsx,"bmaxsy"=>bbmaxsy,"bmaxsz"=>bbmaxsz])

writevtk(trian,"trian",
         celldata=["cellin"=>aggregates,"color"=>colors])

end # module
