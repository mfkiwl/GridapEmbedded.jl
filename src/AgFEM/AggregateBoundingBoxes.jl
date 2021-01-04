function init_bboxes(cell_to_coords)
  # RMK: Assuming first node is min and last node is max of BBox
  [ [cell_to_coords[c][1],cell_to_coords[c][end]] for c in 1:length(cell_to_coords) ]
end

function compute_bboxes!(root_to_agg_bbox,cell_to_root,cell_to_coords)
  for (c,r) in enumerate(cell_to_root)
    ( ( r == 0 ) | ( c == r ) ) && continue
    # RMK: Assuming first node is min and last node is max of BBox
    bbmin = min.(root_to_agg_bbox[r][1].data,cell_to_coords[c][1].data)
    bbmax = max.(root_to_agg_bbox[r][end].data,cell_to_coords[c][end].data)
    root_to_agg_bbox[r] = [bbmin,bbmax]
  end
end

# function extend_bboxes!(root_to_agg_bbox,cell_to_root)
#   for (c,r) in enumerate(cell_to_root)
#     ( ( r == 0 ) | ( c == r ) ) && continue
#     root_to_agg_bbox[c] = root_to_agg_bbox[r]
#   end
# end

function compute_cell_bboxes(model::DiscreteModel,cell_to_root)
  trian = Triangulation(model)
  compute_cell_bboxes(trian,cell_to_root)
end

function compute_cell_bboxes(trian::Triangulation,cell_to_root)
  cell_to_coords = get_cell_coordinates(trian)
  root_to_agg_bbox = init_bboxes(cell_to_coords)
  compute_bboxes!(root_to_agg_bbox,cell_to_root,cell_to_coords)
  # extend_bboxes!(root_to_agg_bbox,cell_to_root)
  root_to_agg_bbox
end

function compute_bbox_dfaces(model::DiscreteModel,cell_to_agg_bbox)
  gt = get_grid_topology(model)
  bboxes = Array{eltype(cell_to_agg_bbox),1}[]
  D = num_dims(gt)
  for d = 1:D-1
    dface_to_Dfaces = get_faces(gt,d,D)
    d_bboxes =
      [ _compute_bbox_dface(dface_to_Dfaces,cell_to_agg_bbox,face) for face in 1:num_faces(gt,d) ]
    bboxes = push!(bboxes,d_bboxes)
  end
  bboxes
end

function _compute_bbox_dface(dface_to_Dfaces,cell_to_agg_bbox,i)
  cells_around_dface_i = getindex(dface_to_Dfaces,i)
  bboxes_around_dface_i = cell_to_agg_bbox[cells_around_dface_i]
  bbmins = [ bboxes_around_dface_i[i][1].data for i in 1:length(bboxes_around_dface_i) ]
  bbmaxs = [ bboxes_around_dface_i[i][2].data for i in 1:length(bboxes_around_dface_i) ]
  [min.(bbmins...),max.(bbmaxs...)]
end


function _compute_cell_to_dface_bboxes(model::DiscreteModel,dbboxes)
  gt = get_grid_topology(model)
  trian = Triangulation(model)
  ctc = get_cell_coordinates(trian)
  bboxes = [ __compute_cell_to_dface_bboxes(gt,ctc,dbboxes,cell) for cell in 1:num_cells(model) ]
  CellPoint(bboxes,trian,PhysicalDomain())
end

function __compute_cell_to_dface_bboxes(gt::GridTopology,ctc,dbboxes,cell::Int)
  cdbboxes = eltype(eltype(eltype(dbboxes)))[]
  D = num_dims(gt)
  Dface_to_0faces = get_faces(gt,D,0)
  for face in getindex(Dface_to_0faces,cell)
    cdbboxes = vcat(cdbboxes,ctc[cell][1],ctc[cell][end])
  end
  for d = 1:D-1
    Dface_to_dfaces = get_faces(gt,D,d)
    for face in getindex(Dface_to_dfaces,cell)
      cdbboxes = vcat(cdbboxes,dbboxes[d][face])
    end
  end
  cdbboxes = vcat(cdbboxes,ctc[cell][1],ctc[cell][end])
end

function compute_cell_to_dface_bboxes(model::DiscreteModel,cell_to_root)
  cbboxes = compute_cell_bboxes(model,cell_to_root)
  dbboxes = compute_bbox_dfaces(model,cbboxes)
  _compute_cell_to_dface_bboxes(model,dbboxes)
end
