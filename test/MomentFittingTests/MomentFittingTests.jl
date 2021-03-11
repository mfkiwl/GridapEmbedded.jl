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

  # function run_test(domain,partition,geom,degree)
  #   bgmodel = CartesianDiscreteModel(domain,partition)
  #   cutgeo = cut(bgmodel,geom)
  #   nodes, moments = compute_cell_moments(cutgeo,degree)
  #   moments = collect(get_array(moments))

  #   # using GridapEmbedded.Interfaces: CUT, IN
  #   # using Gridap.Geometry: get_cell_to_parent_cell

  #   active_model = DiscreteModel(cutgeo,geom,(CUT,IN))
  #   Ωᵃ = Triangulation(active_model)

  #   bgcell_to_inoutcut = compute_bgcell_to_inoutcut(cutgeo,geom)
  #   acell_to_bgcell = get_cell_to_parent_cell(active_model)
  #   acell_to_inoutcut = lazy_map(Reindex(bgcell_to_inoutcut),acell_to_bgcell)
  #   acell_to_point_ptrs = lazy_map(i->(i == CUT ? 1 : 2),acell_to_inoutcut)

  #   # using Gridap.Geometry: get_reffes
  #   # using Gridap.ReferenceFEs: get_polytope
  #   # using Gridap.Integration: Quadrature
  #   quad = map(r->Quadrature(get_polytope(r),2*degree),get_reffes(Ωᵃ))
  #   @assert length(quad) == 1
  #   acell_to_point_vals = [nodes,get_coordinates(quad[1])]

  #   acell_to_weight_vals = moments
  #   push!(acell_to_weight_vals,get_weights(quad[1]))

  #   acell_to_is_cut = findall(lazy_map(i->(i == CUT),acell_to_inoutcut))
  #   num_quads = length(acell_to_weight_vals)
  #   acell_to_weight_ptrs = map(i->(i == IN ? num_quads : 0),acell_to_inoutcut)
  #   acell_to_weight_ptrs[acell_to_is_cut] .= 1:length(acell_to_is_cut)

  #   # using Gridap.Arrays: CompressedArray
  #   acell_to_point = CompressedArray(acell_to_point_vals,acell_to_point_ptrs)
  #   acell_to_weight = CompressedArray(acell_to_weight_vals,acell_to_weight_ptrs)

  #   # using Gridap.Integration: GenericQuadrature
  #   acell_to_quad = [ GenericQuadrature(acell_to_point[i],acell_to_weight[i]) for i in 1:num_quads ]
  #   dΩᵃ = CellQuadrature(acell_to_quad,acell_to_point,acell_to_weight,Ωᵃ,ReferenceDomain())

  #   Ωᶜ = Triangulation(cutgeo)
  #   dΩᶜ = CellQuadrature(Ωᶜ,2*num_dims(bgmodel)*degree)

  #   @test sum(∫(1)*dΩᵃ) - sum(∫(1)*dΩᶜ) + 1 ≈ 1
  #   uᵃ = CellField(x->2*x[1]*x[2],Ωᵃ)
  #   uᶜ = CellField(x->2*x[1]*x[2],Ωᶜ)
  #   @test sum(∫(uᵃ)*dΩᵃ) - sum(∫(uᶜ)*dΩᶜ) + 1 ≈ 1
  # end

  # n = 6
  # degree = 1
  # # partition = (n,n)
  # # domain = (0,1,0,1)
  # # geom = disk(0.42,x0=Point(0.5,0.5))
  # # geom = square(L=0.53,x0=Point(0.5,0.5))
  # partition = (n,n,n)
  # domain = (0,1,0,1,0,1)
  # geom = cube(L=0.53,x0=Point(0.5,0.5,0.5))
  # run_test(domain,partition,geom,degree)

  # # CASE 2D 1
  # n = 6
  # degree = 1
  # partition = (n,n)
  # domain = (0,1,0,1)
  # geom = disk(0.42,x0=Point(0.5,0.5))
  # run_test(domain,partition,geom,degree)

  # # CASE 2D 2
  # n = 6
  # degree = 1
  # partition = (n,n)
  # domain = (0,1,0,1)
  # geom = square(L=0.53,x0=Point(0.5,0.5))
  # run_test(domain,partition,geom,degree)

  # # CASE 3D 1
  # n = 6
  # degree = 1
  # partition = (n,n,n)
  # domain = (0,1,0,1,0,1)
  # geom = sphere(0.42,x0=Point(0.5,0.5,0.5))
  # run_test(domain,partition,geom,degree)

  # # CASE 3D 2
  # n = 6
  # degree = 1
  # partition = (n,n,n)
  # domain = (0,1,0,1,0,1)
  # geom = cube(L=0.53,x0=Point(0.5,0.5,0.5))
  # run_test(domain,partition,geom,degree)

  # CASE 3D 3
  n = 2
  degree = 1
  partition = (n,n,n)
  domain = (0,1,0,1,0,1)
  geom = plane(x0=Point(0.3,0.5,0.5),v=VectorValue(-1.0,0.0,0.0))
  # run_test(domain,partition,geom,degree)

  bgmodel = CartesianDiscreteModel(domain,partition)
  cutgeo = cut(bgmodel,geom)
  nodes, moments = compute_cell_moments(cutgeo,degree)
  moments = collect(get_array(moments))

  # using GridapEmbedded.Interfaces: CUT, IN
  # using Gridap.Geometry: get_cell_to_parent_cell

  active_model = DiscreteModel(cutgeo,geom,(CUT,IN))
  Ωᵃ = Triangulation(active_model)

  bgcell_to_inoutcut = compute_bgcell_to_inoutcut(cutgeo,geom)
  acell_to_bgcell = get_cell_to_parent_cell(active_model)
  acell_to_inoutcut = lazy_map(Reindex(bgcell_to_inoutcut),acell_to_bgcell)
  acell_to_point_ptrs = lazy_map(i->(i == CUT ? 1 : 2),acell_to_inoutcut)

  # using Gridap.Geometry: get_reffes
  # using Gridap.ReferenceFEs: get_polytope
  # using Gridap.Integration: Quadrature
  quad = map(r->Quadrature(get_polytope(r),2*degree),get_reffes(Ωᵃ))
  @assert length(quad) == 1
  acell_to_point_vals = [nodes,get_coordinates(quad[1])]

  acell_to_weight_vals = moments
  push!(acell_to_weight_vals,get_weights(quad[1]))

  acell_to_is_cut = findall(lazy_map(i->(i == CUT),acell_to_inoutcut))
  num_quads = length(acell_to_weight_vals)
  acell_to_weight_ptrs = map(i->(i == IN ? num_quads : 0),acell_to_inoutcut)
  acell_to_weight_ptrs[acell_to_is_cut] .= 1:length(acell_to_is_cut)

  # using Gridap.Arrays: CompressedArray
  acell_to_point = CompressedArray(acell_to_point_vals,acell_to_point_ptrs)
  acell_to_weight = CompressedArray(acell_to_weight_vals,acell_to_weight_ptrs)

  # using Gridap.Integration: GenericQuadrature
  acell_to_quad = [ GenericQuadrature(acell_to_point[i],acell_to_weight[i]) for i in 1:num_quads ]
  dΩᵃ = CellQuadrature(acell_to_quad,acell_to_point,acell_to_weight,Ωᵃ,ReferenceDomain())

  Ωᶜ = Triangulation(cutgeo)
  dΩᶜ = CellQuadrature(Ωᶜ,2*num_dims(bgmodel)*degree)

  @test sum(∫(1)*dΩᵃ) - sum(∫(1)*dΩᶜ) + 1 ≈ 1
  uᵃ = CellField(x->2*x[1]*x[2],Ωᵃ)
  uᶜ = CellField(x->2*x[1]*x[2],Ωᶜ)
  aᵃ = sum(∫(uᵃ)*dΩᵃ)
  aᶜ = sum(∫(uᶜ)*dΩᶜ)
  @test sum(∫(uᵃ)*dΩᵃ) - sum(∫(uᶜ)*dΩᶜ) + 1 ≈ 1

end # module
