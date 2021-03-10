module MomentFittingTests

  using Gridap
  using GridapEmbedded
  using GridapEmbedded.Interfaces: CUT, CUTIN, IN
  using Test

  # CASE 2D 1
  n = 6
  # degree = 2*num_dims(bgmodel)*order
  degree = 1
  partition = (n,n)
  domain = (0,1,0,1)
  geom = disk(0.42,x0=Point(0.5,0.5))
  bgmodel = CartesianDiscreteModel(domain,partition)
  cutgeo = cut(bgmodel,geom)
  Ω = Triangulation(cutgeo)
  vol = sum(get_cell_measure(Ω.a))
  out = sum(sum.(collect(get_array(compute_cell_moments(cutgeo,degree)))))
  @test vol ≈ out

  # CASE 2D 2
  n = 2
  degree = 1
  partition = (n,n)
  domain = (0,1,0,1)
  geom = ! square(L=0.53,x0=Point(0.5,0.5))
  bgmodel = CartesianDiscreteModel(domain,partition)
  cutgeo = cut(bgmodel,geom)
  Ω = Triangulation(cutgeo)
  vol = sum(get_cell_measure(Ω.a))
  out = sum(sum.(collect(get_array(compute_cell_moments(cutgeo,degree)))))
  @test vol ≈ out

  # CASE 3D 1
  n = 6
  degree = 1
  partition = (n,n,n)
  domain = (0,1,0,1,0,1)
  geom = sphere(0.42,x0=Point(0.5,0.5,0.5))
  bgmodel = CartesianDiscreteModel(domain,partition)
  cutgeo = cut(bgmodel,geom)
  Ω = Triangulation(cutgeo)
  vol = sum(get_cell_measure(Ω.a))
  out = sum(sum.(collect(get_array(compute_cell_moments(cutgeo,degree)))))
  @test vol ≈ out

  # CASE 3D 2
  n = 6
  degree = 1
  partition = (n,n,n)
  domain = (0,1,0,1,0,1)
  geom = cube(L=0.53,x0=Point(0.5,0.5,0.5))
  bgmodel = CartesianDiscreteModel(domain,partition)
  cutgeo = cut(bgmodel,geom)
  Ω = Triangulation(cutgeo)
  vol = sum(get_cell_measure(Ω.a))
  out = sum(sum.(collect(get_array(compute_cell_moments(cutgeo,degree)))))
  @test vol ≈ out

  # CASE 3D 3
  n = 2
  degree = 1
  partition = (n,n,n)
  domain = (0,1,0,1,0,1)
  geom = plane(x0=Point(0.3,0.5,0.5),v=VectorValue(-1.0,0.0,0.0))
  bgmodel = CartesianDiscreteModel(domain,partition)
  cutgeo = cut(bgmodel,geom)
  Ω = Triangulation(cutgeo)
  vol = sum(get_cell_measure(Ω.a))
  out = sum(sum.(collect(get_array(compute_cell_moments(cutgeo,degree)))))
  @test vol ≈ out

  # Ω = Triangulation(cutgeo)
  # writevtk(Ω,"trian")
  # Γᵉ = EmbeddedBoundary(cutgeo)
  # writevtk(Γᵉ,"skeletone")
  # Λ  = GhostSkeleton(cutgeo)
  # cutgeo_facets = cut_facets(cutgeo.bgmodel,cutgeo.geo)
  # Γᶠ = SkeletonTriangulation(cutgeo_facets,Λ,cutgeo.geo,CUTIN)
  # writevtk(Γᶠ,"skeletonf")
  # Γᵇ = SkeletonTriangulation(cutgeo_facets,Λ,cutgeo.geo,IN)
  # if num_cells(Γᵇ) > 0
  #   writevtk(Γᵇ,"skeletonb")
  # end
  # Γᵒ = BoundaryTriangulation(cutgeo_facets,CUTIN)
  # if num_cells(Γᵒ) > 0
  #   writevtk(Γᵒ,"skeletono")
  # end

end # module
