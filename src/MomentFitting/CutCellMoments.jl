function compute_cell_moments(model::RestrictedDiscreteModel,
                              facet_moments::DomainContribution)

  cell_to_parent_cell = get_cell_to_parent_cell(model)
  fi = [ testitem(array) for (trian,array) in facet_moments.dict ]
  li = map(length,fi)
  @assert all(li .== first(li))

  cut_cell_to_moments = Fill(zero(first(fi)),length(cell_to_parent_cell))

  for (trian,array) in facet_moments
    # compute_cell_moments(trian,array) # to implement
  end

  cut_cell_to_moments

end
