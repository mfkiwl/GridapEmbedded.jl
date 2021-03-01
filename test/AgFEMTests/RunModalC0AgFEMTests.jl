using MPIClusterManagers, Distributed
import MPI

manager = MPIClusterManagers.start_main_loop(MPI_TRANSPORT_ALL)

@everywhere function pidandhost()
  host, pid = (gethostname(), getpid())
  println("Hello from process $(pid) on host $(host)!")
end

@everywhere include("test/AgFEMTests/ModalC0AgFEMTests.jl")
@everywhere using .ModalC0AgFEMTests
@everywhere using DrWatson
@everywhere @quickactivate "GridapEmbedded"
@everywhere using GridapPETSc
@everywhere GridapPETSc.Init(["-ksp_type", "cg",
                              "-ksp_rtol", "$tol",
                              "-ksp_max_it", "$maxits",
                              "-ksp_norm_type", "unpreconditioned",
                              "-pc_type","gamg",
                              "-pc_gamg_type","agg",
                              "-pc_gamg_esteig_ksp_type","cg",
                              "-mg_levels_esteig_ksp_type","cg",
                              "-mg_coarse_sub_pc_type","cholesky",
                              "-mg_coarse_sub_pc_factor_mat_ordering_type","nd",
                              "-pc_gamg_process_eq_limit","50",
                              "-pc_gamg_square_graph","0",
                              "-pc_gamg_agg_nsmooths","1"])

params = Dict(
  :n => [6,12,24,48,96],
  :k => [1,2,3,4,5],
  :d => [2],
  :t => [0,1],
  :s => [0,1],
  :g => [0,1]
)
dicts = dict_list(params)

@everywhere @info "Training"
@everywhere compute(6,1,2,1,0,0)
@everywhere compute(6,3,2,1,0,0)
@everywhere @info "Producing"
pmap(compute_and_save,dicts)

@everywhere GridapPETSc.Finalize()

MPIClusterManagers.stop_main_loop(manager)
