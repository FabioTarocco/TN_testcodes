
#----------------MAIN------------------
using ITensors
using PlotlyJS

using LinearAlgebra
include("mutual_info.jl")
include("rdm.jl")
function job()


    Ns = 4;
    N = Int64(.5 * Ns * (Ns+1));

    diag = true;

    sites = siteinds("S=1/2", N, conserve_qns = false);
    dlattice = [(0, 4), (1, 4), (0, 1), (1, 5), (2, 5), (1, 2), (2, 6), (3, 6), (2, 3), (4, 7), (5, 7), (4, 5), (5, 8), (6, 8), (5, 6), (7, 9), (8, 9), (7, 8)];




    for b in dlattice
        println(b)
    end

    triangle_sumOp = OpSum()

    for b in dlattice
        triangle_sumOp +=  "Sz", b[1]+1, "Sz", b[2]+1
        triangle_sumOp +=  "Sx", b[1]+1, "Sx", b[2]+1
        triangle_sumOp +=  "Sy", b[1]+1, "Sy", b[2]+1
    end

    triangle_H = MPO(triangle_sumOp, sites);
    @show triangle_H;
    @show triangle_sumOp


    
    if diag == true
        prod_H = prod(triangle_H);
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
        @show (eig_vals[2])
        @show (eig_vals[3])
        #Initial configuration of the sites of the system
        #initial_state = [isodd(n) ? "Up" : "Dn" for n=1:N]
    end
    initial_state = [isodd(n) ? "Up" : "Dn" for n=1:N]
    psi_init_8 = randomMPS(sites, initial_state);


    sweeps = Sweeps(10)
    setmaxdim!(sweeps, 8,16,50,100,200,400,1000,2000)
    setcutoff!(sweeps, 1e-10)
    setnoise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0)
    energy_DMRG, psi_DMRG = dmrg(triangle_H, psi_init_8, sweeps);
    @show energy_DMRG

    if maxlinkdim(psi_DMRG) > 8
        psi_DMRG_trunc = psi_DMRG
        truncate!(psi_DMRG_trunc, maxlinkdim=8)
    end

    @show energy_DMRG - eig_vals[1]
    @show inner(psi_DMRG_trunc, triangle_H, psi_DMRG_trunc) - eig_vals[1]

    full_DM = reduced_rho_matrix(psi_DMRG);
    I_matrix = MI_diag(full_DM);
    maxC = findmax(I_matrix);
    I_matrixN = I_matrix ./ maxC[1];

    plot(heatmap(z=I_matrixN, colorscale = "Viridis"))
    bond_dim = maxlinkdim(psi_DMRG);
    path_N= string("heisenberg_triangle_$(Ns)_chi_$(bond_dim).png");
    savefig(plot(heatmap(z=I_matrixN, colorscale = "Viridis")), path_N)
end