function init_bboxes(cell_to_coords)
  T = eltype(eltype(cell_to_coords))
  # RMK: Assuming first node is min and last node is max of BBox
  [ (cell_to_coords[c][1],cell_to_coords[c][end]) for c in 1:length(cell_to_coords) ]
end

function compute_bboxes!(root_to_agg_bbox,cell_to_root,cell_to_coords)
  for (c,r) in enumerate(cell_to_root)
    ( ( r == 0 ) | ( c == r ) ) && continue
    # RMK: Assuming first node is min and last node is max of BBox
    bbmin = min.(root_to_agg_bbox[r][1].data,cell_to_coords[c][1].data)
    bbmax = max.(root_to_agg_bbox[r][end].data,cell_to_coords[c][end].data)
    root_to_agg_bbox[r] = (bbmin,bbmax)
  end
end

# function extend_bboxes!(root_to_agg_bbox,cell_to_root)
#   for (c,r) in enumerate(cell_to_root)
#     ( ( r == 0 ) | ( c == r ) ) && continue
#     root_to_agg_bbox[c] = root_to_agg_bbox[r]
#   end
# end

function compute_aggregate_bboxes(model::DiscreteModel,cell_to_root)
  trian = Triangulation(model)
  compute_aggregate_bboxes(trian,cell_to_root)
end

function compute_aggregate_bboxes(trian::Triangulation,cell_to_root)
  cell_to_coords = collect(get_cell_coordinates(trian))
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
  (min.(bbmins...),max.(bbmaxs...))
end

function compute_anchor_nfaces(model::DiscreteModel)
  gt = get_grid_topology(model)
  vc = get_vertex_coordinates(gt)
  anchors = Array{eltype(vc),1}[]
  for d = 1:num_dims(gt)-1
    dface_to_0faces = get_faces(gt,d,0)
    d_anchors =
      [ vc[dface_to_0faces.data[dface_to_0faces.ptrs[i]]] for i in 1:num_faces(gt,d) ]
    anchors = push!(anchors,d_anchors)
  end
  anchors
end
