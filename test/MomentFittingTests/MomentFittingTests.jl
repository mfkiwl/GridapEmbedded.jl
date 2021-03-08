module MomentFittingTests

  using Gridap
  using Gridap.CellData
  using Gridap.CellData: SkeletonPair
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
  Cᵉ = compute_hyperplane_coeffs(Γᵉ)
  writevtk(Γᵉ,"skeletone",celldata=["Cᵉ"=>Cᵉ])

  Λ  = GhostSkeleton(cutgeo)
  Γᶠ = SkeletonTriangulation(cutgeo_facets,Λ,geom,CUTIN)
  C⁺ = compute_hyperplane_coeffs(Γᶠ,:⁺)
  C⁻ = compute_hyperplane_coeffs(Γᶠ,:⁻)
  writevtk(Γᶠ.⁺,"skeletonf")
  # writevtk(Γᶠ.⁺,"skeletonf",celldata=["plus"=>C⁺]) # Cannot draw field

  bgfacet_to_inoutcut = compute_bgfacet_to_inoutcut(cfgeom,geom)
  bgfacet_to_mask = lazy_map( a -> a == IN, bgfacet_to_inoutcut )
  Γᵇ = BoundaryTriangulation(cutbgm,bgfacet_to_mask)
  Cᵇ = compute_hyperplane_coeffs(Γᵇ)
  writevtk(Γᵇ,"skeletonb")
  # writevtk(Γᵇ,"skeletonb",celldata=["coeffs"=>Cᵇ]) # Cannot draw field

  dΓᵉ = Measure(Γᵉ,2*D*order)
  dΓᶠ = Measure(Γᶠ,2*D*order)
  dΓᵇ = Measure(Γᵇ,2*D*order)

  T = Float64
  b = MonomialBasis{2}(T,order)
  i = Matrix{eltype(T)}(I,length(b),length(b))
  l = linear_combination(i,b)

  m = Fill(l,num_cells(Ωᵇ))
  v = GenericCellField(m,Ωᵇ,ReferenceDomain())

  cᵉ = CellField(Cᵉ,Γᵉ)
  c⁺ = GenericCellField(C⁺,Γᶠ,ReferenceDomain())
  # c⁻ = GenericCellField(C⁻,Γᶠ,ReferenceDomain())
  cᵇ = GenericCellField(Cᵇ,Γᵇ,ReferenceDomain())

  # # Compute integrals on the boundary facets of the polytope
  # Iᵉ = ∫( v   )*dΓᵉ # Embedded facets
  # I⁺ = ∫( v.⁺ )*dΓᶠ # Cut interior (fitted) facets (plus)
  # # I⁻ = ∫( v.⁻ )*dΓᶠ # Cut interior (fitted) facets (minus)
  # Iᵇ = ∫( v   )*dΓᵇ # Interior (fitted) facets

  # Jᵉ = ∫( cᵉ*v   )*dΓᵉ # Embedded facets
  # J⁺ = ∫( c⁺*v.⁺ )*dΓᶠ # Cut interior (fitted) facets (plus)
  # # J⁻ = ∫( c⁺*v.⁻ )*dΓᶠ # Cut interior (fitted) facets (minus)
  # Jᵇ = ∫( cᵇ*v   )*dΓᵇ # Interior (fitted) facets

  J = ∫(cᵉ*v)*dΓᵉ + ∫(c⁺*v.⁺)*dΓᶠ + ∫(cᵇ*v)*dΓᵇ
  out = compute_cell_moments(cutbgm,J)

end # module
