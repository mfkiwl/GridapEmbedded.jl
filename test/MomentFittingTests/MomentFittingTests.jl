module MomentFittingTests

  using Gridap
  using GridapEmbedded

  const R = 2.52
  n = 6
  order = 1

  # geom = square(L=0.63,x0=Point(0.5,0.5))
  geom = disk(R,x0=Point(3.0,3.0))
  partition = (n,n)
  domain = (0,6,0,6)

  bgmodel = CartesianDiscreteModel(domain,partition)

  cutgeo = cut(bgmodel,geom)
  degree = 1 # 2*order

  out = sum.(collect(compute_cell_moments(cutgeo,degree)))

end # module
