using ITensors
using PlotlyJS

using LinearAlgebra
include("mutual_info.jl")
include("rdm.jl")

function map()
    color = 2;
    sites = 6;
    total_sites = color*sites;
    hopping = [(0,1),(1,2),(2,3),(3,4),(4,5)];
    t_term = 1.0;
    U_term = 8.0 ;
    mu_term = 4.0;
    diag = false

    
    s = siteinds("Fermion",total_sites);
    hubbard_sumOp = OpSum();
    pairs = [];
    for p in hopping
        push!(pairs, (p[1]+1, p[2]+1))
    end

    for p in pairs
        hubbard_sumOp -= t_term, "cdag", p[1]*2-1, "c", p[2]*2-1;
        hubbard_sumOp -= t_term, "cdag", p[2]*2-1, "c", p[1]*2-1;
        hubbard_sumOp -= t_term, "cdag", p[1]*2, "c", p[2]*2;
        hubbard_sumOp -= t_term, "cdag", p[2]*2, "c", p[1]*2;
    end
    for i in 1:sites
        hubbard_sumOp -= (2*mu_term + U_term), "cdag", i*2-1, "cdag", i*2, "c", i*2-1, "c", i*2
    end
    for i in 1:sites
        hubbard_sumOp -= mu_term, "cdag", i*2-1, "c", i*2-1
        hubbard_sumOp -= mu_term, "cdag", i*2, "c", i*2
    end

    const_energy = 4*mu_term*sites;

    hubbard_MPO = MPO(hubbard_sumOp, s);
    @show hubbard_sumOp;
    @show hubbard_MPO;

    if diag == true
        prod_H = prod(hubbard_MPO);
        @show inds(prod_H);

        list_inds_prime = []
        for n in 1:2*total_sites
            if isodd(n)
                push!(list_inds_prime, inds(prod_H)[n])
            end
        end;
        @show list_inds_prime;

        list_inds_noprime = []
        for n in 1:2*total_sites
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
    end
    fermion_state = [isodd(n) ? 1 : 2 for n=1:total_sites];
    init_half_filling = randomMPS(s, fermion_state);

    sweeps = Sweeps(10);
    setmaxdim!(sweeps, 10,20,50,100,200,400,1000);
    setnoise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0);
    energy_DMRG, psi_DMRG = dmrg(hubbard_MPO, init_half_filling, sweeps);
    @show energy_DMRG + const_energy;
    # ED = -24.92191748244516 

end