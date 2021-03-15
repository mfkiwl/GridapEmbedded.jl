module MomentFittingTests

  using Gridap
  using GridapEmbedded
  using Test

  using GridapEmbedded.Interfaces: CUT, IN

  using Gridap.Geometry: get_cell_to_parent_cell

  using Gridap.Geometry: get_reffes
  using Gridap.ReferenceFEs: get_polytope
  using Gridap.Integration: Quadrature

  using Gridap.Arrays: CompressedArray

  using Gridap.Integration: GenericQuadrature

  function run_test(domain,partition,geom,degree)
    bgmodel = CartesianDiscreteModel(domain,partition)
    cutgeo = cut(bgmodel,geom)
    active_model = DiscreteModel(cutgeo,geom,(CUT,IN))
    Ωᵃ = Triangulation(active_model)
    dΩᵃ = Measure(MomentFittingQuad(Ωᵃ,cutgeo,degree))
    Ωᶜ = Triangulation(cutgeo)
    dΩᶜ = Measure(Ωᶜ,2*num_dims(bgmodel)*degree)
    @test sum(∫(1)*dΩᵃ) - sum(∫(1)*dΩᶜ) + 1 ≈ 1
    uᵃ = CellField(x->2*x[1]*x[2],Ωᵃ)
    uᶜ = CellField(x->2*x[1]*x[2],Ωᶜ)
    @test sum(∫(uᵃ)*dΩᵃ) - sum(∫(uᶜ)*dΩᶜ) + 1 ≈ 1
  end

  # CASE 2D 1
  n = 6
  degree = 1
  partition = (n,n)
  domain = (0,1,0,1)
  geom = disk(0.42,x0=Point(0.5,0.5))
  run_test(domain,partition,geom,degree)

  # CASE 2D 2
  n = 6
  degree = 1
  partition = (n,n)
  domain = (0,1,0,1)
  geom = square(L=0.53,x0=Point(0.5,0.5))
  run_test(domain,partition,geom,degree)

  # CASE 3D 1
  n = 6
  degree = 1
  partition = (n,n,n)
  domain = (0,1,0,1,0,1)
  geom = sphere(0.42,x0=Point(0.5,0.5,0.5))
  run_test(domain,partition,geom,degree)

  # CASE 3D 2
  n = 6
  degree = 1
  partition = (n,n,n)
  domain = (0,1,0,1,0,1)
  geom = cube(L=0.53,x0=Point(0.5,0.5,0.5))
  run_test(domain,partition,geom,degree)

  # CASE 3D 3
  n = 2
  degree = 1
  partition = (n,n,n)
  domain = (0,1,0,1,0,1)
  geom = plane(x0=Point(0.3,0.5,0.5),v=VectorValue(-1.0,0.0,0.0))
  run_test(domain,partition,geom,degree)

end # module
