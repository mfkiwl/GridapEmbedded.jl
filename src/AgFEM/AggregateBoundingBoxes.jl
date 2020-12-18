function init_bboxes(cell_to_coords,cell_to_root)
  T = eltype(eltype(cell_to_coords))
  # RMK: Assuming first node is min and last node is max of BBox
  root_bbox(c,r,T) =
    c == r ? (cell_to_coords[c][1],cell_to_coords[c][end]) : (zero(T),zero(T))
  [ root_bbox(c,r,T) for (c,r) in enumerate(cell_to_root) ]
end

function update_bboxes!(root_to_agg_bbox,cell_to_coords,cell_to_root)
  for (c,r) in enumerate(cell_to_root)
    ( ( r == 0 ) | ( c == r ) ) && continue
    # RMK: Assuming first node is min and last node is max of BBox
    bbmin = min.(root_to_agg_bbox[r][1].data,cell_to_coords[c][1].data)
    bbmax = max.(root_to_agg_bbox[r][end].data,cell_to_coords[c][end].data)
    root_to_agg_bbox[r] = (bbmin,bbmax)
  end
end

function compute_aggregate_bboxes(trian,cell_to_root)
  cell_to_coords = collect(get_cell_coordinates(trian))
  root_to_agg_bbox = init_bboxes(cell_to_coords,cell_to_root)
  update_bboxes!(root_to_agg_bbox,cell_to_coords,cell_to_root)
  root_to_agg_bbox
end
