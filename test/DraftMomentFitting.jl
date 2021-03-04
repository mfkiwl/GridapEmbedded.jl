module DraftMomentFitting

  using Gridap
  using Gridap.CellData
  using Gridap.Polynomials
  using Gridap.Fields: linear_combination
  using GridapEmbedded
  using GridapEmbedded.Interfaces: CUT, CUTIN, IN

  using FillArrays
  import LinearAlgebra: I

  const R = 2.52

  n = 6
  order = 1

  # geom = square(L=0.63,x0=Point(0.5,0.5))
  geom = disk(R,x0=Point(3.0,3.0))
  partition = (n,n)
  domain = (0,6,0,6)

  bgmodel = CartesianDiscreteModel(domain,partition)
  D = num_dims(bgmodel)
  h = (domain[2]-domain[1])/n

  cutgeo = cut(bgmodel,geom)
  cutgeo_facets = cut_facets(bgmodel,geom)

  cutbgm = DiscreteModel(cutgeo,geom,CUT)
  cgeom  = cut(cutbgm,geom)
  cfgeom = cut_facets(cutbgm,geom)

  Ωᵇ = Triangulation(bgmodel)
  Γᵉ = EmbeddedBoundary(cutgeo)
  # writevtk(Γᵉ,"skeletone")

  Λ  = GhostSkeleton(cutgeo)
  Γᶠ = SkeletonTriangulation(cutgeo_facets,Λ,geom,CUTIN)
  # writevtk(Γᶠ,"skeletonf")

  bgfacet_to_inoutcut = compute_bgfacet_to_inoutcut(cfgeom,geom)
  bgfacet_to_mask = lazy_map( a -> a == IN, bgfacet_to_inoutcut )
  Γᵇ = BoundaryTriangulation(cutbgm,bgfacet_to_mask)
  # writevtk(Γᵇ,"skeletonb")

  dΓᵉ = CellQuadrature(Γᵉ,2*D*order)
  dΓᶠ = CellQuadrature(Γᶠ,2*D*order)
  dΓᵇ = CellQuadrature(Γᵇ,2*D*order)

  T = Float64
  b = MonomialBasis{2}(T,order)
  i = Matrix{eltype(T)}(I,length(b),length(b))
  l = linear_combination(i,b)

  m = Fill(l,num_cells(Ωᵇ))
  v = GenericCellField(m,Ωᵇ,ReferenceDomain())

  # Compute integrals on the boundary facets of the polytope
  Iᵉ = ∫( v   )*dΓᵉ # Embedded facets
  Iᶠ = ∫( v.⁺ )*dΓᶠ # Cut interior (fitted) facets
  Iᵇ = ∫( v   )*dΓᵇ # Interior (fitted) facets

end # module
