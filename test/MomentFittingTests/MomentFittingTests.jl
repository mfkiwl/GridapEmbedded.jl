module MomentFittingTests

  using Gridap
  using GridapEmbedded

  using GridapEmbedded.Interfaces: CUTIN, IN

  n = 2
  order = 1

  geom = ! square(L=0.53,x0=Point(0.5,0.5))
  # geom = disk(0.42,x0=Point(0.5,0.5))
  partition = (n,n)
  domain = (0,1,0,1)

  # geom = cube(L=0.53,x0=Point(0.5,0.5,0.5))
  # geom = sphere(0.42,x0=Point(0.5,0.5,0.5))
  # partition = (n,n,n)
  # domain = (0,1,0,1,0,1)

  bgmodel = CartesianDiscreteModel(domain,partition)

  cutgeo = cut(bgmodel,geom)
  cutgeo_facets = cut_facets(bgmodel,geom)

  Ω = Triangulation(cutgeo)
  writevtk(Ω,"trian")
  Γᵉ = EmbeddedBoundary(cutgeo)
  writevtk(Γᵉ,"skeletone")
  Λ  = GhostSkeleton(cutgeo)
  cutgeo_facets = cut_facets(cutgeo.bgmodel,cutgeo.geo)
  Γᶠ = SkeletonTriangulation(cutgeo_facets,Λ,cutgeo.geo,CUTIN)
  writevtk(Γᶠ,"skeletonf")
  Γᵇ = SkeletonTriangulation(cutgeo_facets,Λ,cutgeo.geo,IN)
  if num_cells(Γᵇ) > 0
    writevtk(Γᵇ,"skeletonb")
  end
  bgfacet_to_inoutcut = compute_bgfacet_to_inoutcut(cutgeo_facets,geom)
  bgfacet_to_mask = lazy_map( a -> a == IN, bgfacet_to_inoutcut )
  Γᵒ = BoundaryTriangulation(bgmodel,bgfacet_to_mask)
  if num_cells(Γᵒ) > 0
    writevtk(Γᵒ,"skeletono")
  end

  # degree = 2*num_dims(bgmodel)*order
  degree = 1

  out = sum.(collect(compute_cell_moments(cutgeo,degree)))

end # module
