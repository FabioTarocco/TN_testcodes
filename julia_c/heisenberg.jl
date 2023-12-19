#----------------MAIN------------------
using ITensors
using PlotlyJS

using LinearAlgebra
include("mutual_info.jl")
include("rdm.jl")
function job()
    A = [
        0 1 0 0 0 1 0 0 0;
        0 0 1 0 1 0 0 0 0;
        0 0 0 1 0 0 0 0 0;
        0 0 0 0 1 0 0 0 1;
        0 0 0 0 0 1 0 1 0;
        0 0 0 0 0 0 1 0 0;
        0 0 0 0 0 0 0 1 0;
        0 0 0 0 0 0 0 0 1;
        0 0 0 0 0 0 0 0 0;
    ]

    num_nodes = size(A, 1)

    edge_pairs = []

    for i in 1:num_nodes
        for j in 1:num_nodes
            if A[i, j] == 1
                push!(edge_pairs, (i, j))
            end
        end
    end

    for pair in edge_pairs
        println(pair)
    end



    Nx = 3;
    Ny = 4;
    N = Nx * Ny;

    mps_mode = 0 #0=std mps, 1=snake
    lattice_type = 0 # 0=open boundary, 1=periodic on Y
    heisenberg_mode = 0 #0=ZZXXYY, 1=+--+ZZ
    diag = true

    sites = siteinds("S=1/2", N, conserve_qns = false)
    #sites = siteind("Qubit", N, conserve_qns = false)

    if lattice_type == 0
        dlattice = square_lattice(Nx, Ny; yperiodic=false)
    elseif lattice_type == 1
        dlattice = square_lattice(Nx, Ny; yperiodic=true)
    end

    for b in dlattice
        println(b)
    end

    heisenberg_2D_sumOp = OpSum();

    if mps_mode == 1
        if heisenberg_mode == 0 

            for pair in edge_pairs
                heisenberg_2D_sumOp +=  "Sz", pair[1], "Sz", pair[2]
                heisenberg_2D_sumOp +=  "Sx", pair[1], "Sx", pair[2]
                heisenberg_2D_sumOp +=  "Sy", pair[1], "Sy", pair[2]
            end
        elseif heisenberg_mode == 1
            for pair in edge_pairs
                heisenberg_2D_sumOp +=       "Sz", pair[1], "Sz", pair[2]
                heisenberg_2D_sumOp += 1/2,  "S+", pair[1], "S-", pair[2]
                heisenberg_2D_sumOp += 1/2,  "S-", pair[1], "S+", pair[2]
            end
        end
    elseif mps_mode == 0

        if heisenberg_mode == 0 
            for b in dlattice
                heisenberg_2D_sumOp +=  "Sz", b.s1, "Sz", b.s2
                heisenberg_2D_sumOp +=  "Sx", b.s1, "Sx", b.s2
                heisenberg_2D_sumOp +=  "Sy", b.s1, "Sy", b.s2
            end
        elseif heisenberg_mode == 1
            for b in dlattice
                heisenberg_2D_sumOp +=       "Sz", b.s1, "Sz", b.s2
                heisenberg_2D_sumOp += 1/2,  "S+", b.s1, "S-", b.s2
                heisenberg_2D_sumOp += 1/2,  "S-", b.s1, "S+", b.s2
            end
        end
    end
    #MPO of the 2D Heisenberg Hamiltonian

    heisenberg_2D_H = MPO(heisenberg_2D_sumOp, sites);
    @show heisenberg_2D_H;
    @show heisenberg_2D_sumOp

    if diag == true
        prod_H = prod(heisenberg_2D_H);
        @show inds(prod_H);

        list_inds_prime = []
        for n in 1:2*N
            if isodd(n)
                push!(list_inds_prime, inds(prod_H)[n])
            end
        end;
        @show list_inds_prime;

        list_inds_noprime = []
        for n in 1:2*N
            if iseven(n)
                push!(list_inds_noprime, inds(prod_H)[n])
            end
        end;
        @show list_inds_noprime;
        combiner_prime = combiner(list_inds_prime);
        combiner_noprime = combiner(list_inds_noprime);

        matrix_H = prime(combiner_prime, "Link")*prod_H*combiner_noprime;
        @show inds(matrix_H);
        @show size(matrix_H);

        eig_vals, eig_vec = eigen(matrix_H, inds(matrix_H)[1], inds(matrix_H)[2]);
        @show (eig_vals[1])
        #Initial configuration of the sites of the system
        #initial_state = [isodd(n) ? "Up" : "Dn" for n=1:N]
    end
    @show real(eig_vals[1])

    initial_state = [isodd(n) ? "Up" : "Dn" for n=1:N]
    psi_init_8 = randomMPS(sites, initial_state);


    sweeps = Sweeps(10)
    setmaxdim!(sweeps, 16,16,16,16,16,16,16,16,16,16)
    setcutoff!(sweeps, 1e-10)
    setnoise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0)
    energy_DMRG, psi_DMRG = dmrg(heisenberg_2D_H, psi_init_8, sweeps);
    @show energy_DMRG

    @show maxlinkdim(psi_DMRG)
    if maxlinkdim(psi_DMRG) > 8
        println("ENTRA")
        psi_DMRG_trunc = truncate(maxlinkdim=16)
    end

    @show energy_DMRG - eig_vals[1]
    @show inner(psi_DMRG, heisenberg_2D_H, psi_DMRG) - eig_vals[1]

    #=
    full_DM_trunc = reduced_rho_matrix(psi_DMRG_trunc);
    I_matrix_trunc = MI(full_DM_trunc);
    maxC_trunc = findmax(I_matrix_trunc);
    I_matrixN_trunc = I_matrix_trunc ./ maxC_trunc[1];

    plot(heatmap(z=I_matrixN_trunc, colorscale = "Viridis"))
    path_T= string("heisenberg%dx%d_truncaded_chi_%d_test.png",Nx,Ny,maxlinkdim(psi_DMRG))
    savefig(plot(heatmap(z=I_matrixN_trunc, colorscale = "Viridis")), path_T)
    =#
    full_DM = reduced_rho_matrix(psi_DMRG);
    I_matrix = MI_diag(full_DM);
    maxC = findmax(I_matrix);
    I_matrixN = I_matrix ./ maxC[1];

    plot(heatmap(z=I_matrixN, colorscale = "Viridis"))
    bond_dim = maxlinkdim(psi_DMRG);
    path_N= string("heisenberg$(Nx)x$(Ny)_chi_$(bond_dim)_truncated.png");
    savefig(plot(heatmap(z=I_matrixN, colorscale = "Viridis")), path_N)


    I_diag_matrix = MI_diag(full_DM);
    I_diag_matrix = I_diag_matrix ./ findmax(I_diag_matrix)[1];
    plot(heatmap(z=I_diag_matrix, colorscale = "Viridis"))



    #=
    mi_plot = plot_MI_coupling(I_matrixN, size(I_matrixN)[1])

    @show mi_plot

    mi_plot
    p = plot()
    i = 0
    for k ∈ keys(mi_plot)
        scatter!(p, mi_plot[k],i)
    end
    =#




end
