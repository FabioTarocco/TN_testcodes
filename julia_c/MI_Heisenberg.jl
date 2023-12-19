#----------------MAIN------------------
using ITensors
using PlotlyJS

using LinearAlgebra
include("mutual_info.jl")
include("rdm.jl")
global save = true
function diag(op, N)

    prod_H = prod(op);
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
    return eig_vals, eig_vec
end

function map_squared(nx, ny, ED_ref, trunc_th, is_trunc, delta)

    Nx = nx;
    Ny = ny;
    N = Nx * Ny;

    heisenberg_mode = 0;

    dlattice = square_lattice(Nx, Ny; yperiodic=false)

    sites = siteinds("S=1/2", N, conserve_qns = false)

    heisenberg_2D_sumOp = OpSum();

    if heisenberg_mode == 0 
        for b in dlattice
            heisenberg_2D_sumOp +=  "Sz", b.s1, "Sz", b.s2
            heisenberg_2D_sumOp +=  "Sx", b.s1, "Sx", b.s2
            heisenberg_2D_sumOp +=  delta,"Sy", b.s1, "Sy", b.s2
        end
    elseif heisenberg_mode == 1
        for b in dlattice
            heisenberg_2D_sumOp +=       "Sz", b.s1, "Sz", b.s2
            heisenberg_2D_sumOp += 1/2,  "S+", b.s1, "S-", b.s2
            heisenberg_2D_sumOp += 1/2,  "S-", b.s1, "S+", b.s2
        end
    end

    #MPO of the 2D Heisenberg Hamiltonian

    heisenberg_2D_H = MPO(heisenberg_2D_sumOp, sites);
    @show heisenberg_2D_H;
    @show heisenberg_2D_sumOp

    if ED_ref
        eig_vals, eig_vec = diag(heisenberg_2D_H, N)
        @show eig_vals[1]
    end

    initial_state = [isodd(n) ? "Up" : "Dn" for n=1:N]
    init_MPS = randomMPS(sites, initial_state);


    sweeps = Sweeps(10)
    if is_trunc
        setmaxdim!(sweeps, trunc_th)
    else
        setmaxdim!(sweeps, 10,20,50,100,200,500,800,1000)
    end

    setcutoff!(sweeps, 1e-10)
    setnoise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0)
    energy_DMRG, psi_DMRG = dmrg(heisenberg_2D_H, init_MPS, sweeps);
    @show energy_DMRG

    if ED_ref
        @show energy_DMRG - eig_vals[1]
        delta = inner(psi_DMRG', heisenberg_2D_H, psi_DMRG) - eig_vals[1]
    end

    full_DM = reduced_rho_matrix(psi_DMRG);
    I_matrix,norm_factor, mi_pairs = MI_diag(full_DM);

    plot(heatmap(z=I_matrix, colorscale = "Viridis"))
    bond_dim = maxlinkdim(psi_DMRG);
    if save
        if delta == 1.0
            path_N= string("../Image_results/MI_Heiseberg_$(Nx)x$(Ny)_chi_$(bond_dim).png");
        else
            path_N= string("../Image_results/MI_Heiseberg_$(Nx)x$(Ny)_dy0_chi_$(bond_dim).png")
        end
        savefig(plot(heatmap(z=I_matrix, colorscale = "Viridis")), path_N)
    end
    
    if is_trunc
        println("Truncation threshold: $(trunc_th)")
    end
    println("SYSTEM:\n$(Nx)x$(Ny)\n")
    println("Julia lattice indices")
    for b in dlattice
        println(b)
    end 
    println("Python lattice indices")
    for b in dlattice
        println("Coupling:",(b.s1 -1, b.s2 -1))
    end
    if ED_ref
        println("Exact Diagonalization Energy: $(eig_vals[1])")
    end
    println("DMRG energy: $(energy_DMRG)")
    println("Normalization factor: $(norm_factor)")
    println(mi_pairs)
end



function map_triangle(ED_ref, trunc_th, is_trunc)


    Ns = 4;
    N = Int64(.5 * Ns * (Ns+1));


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

    triangle_MPO = MPO(triangle_sumOp, sites);
    @show triangle_MPO;
    @show triangle_sumOp
    @show size(triangle_sumOp)[1] == size(dlattice)[1]*3;


    if ED_ref
        eig_vals, eig_vec = diag(triangle_MPO, N)
        @show eig_vals[1]
    end

    initial_state = [isodd(n) ? "Up" : "Dn" for n=1:N]
    init_MPS = randomMPS(sites, initial_state);


    sweeps = Sweeps(10)
    if is_trunc
        setmaxdim!(sweeps, trunc_th)
    else
        setmaxdim!(sweeps, 10,20,50,100,200,500,800,1000)
    end

    setcutoff!(sweeps, 1e-10)
    setnoise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0)

    energy_DMRG, psi_DMRG = dmrg(triangle_MPO, init_MPS, sweeps);
    @show energy_DMRG

    if ED_ref
        @show energy_DMRG - eig_vals[1]
        delta = inner(psi_DMRG', triangle_MPO, psi_DMRG) - eig_vals[1]
    end

    full_DM = reduced_rho_matrix(psi_DMRG);
    I_matrix,norm_factor, mi_pairs = MI_diag(full_DM);

    plot(heatmap(z=I_matrix, colorscale = "Viridis"))
    bond_dim = maxlinkdim(psi_DMRG);

    if save
    path_N= string("../Image_results/MI_Heiseberg_Triangle_$(N)_chi_$(bond_dim).png");
    savefig(plot(heatmap(z=I_matrix, colorscale = "Viridis")), path_N)
    end

    
    if is_trunc
        println("Truncation threshold: $(trunc_th)")
    end
    println("SYSTEM:\nHeisenberg triangle lattice $(N) sites\n")
    println("Julia lattice indices")
    for b in dlattice
        println("Coupling:",(b[1]+1,b[2]+1))
    end 
    println("Python lattice indices")
    for b in dlattice
        println("Coupling:",(b[1], b[2]))
    end
    if ED_ref
        println("Exact Diagonalization Energy: $(eig_vals[1])")
    end
    println("DMRG energy: $(energy_DMRG)")
    println("Normalization factor: $(norm_factor)")
    println(mi_pairs)
end



function map_hex(ED_ref, trunc_th, is_trunc)
    N = 9;

    dlattice = [(0, 3), (1, 3), (0, 1), (1, 4), (2, 4), (1, 2), (2, 5), (3, 6), (4, 6), (3, 4), (6, 7), (4, 7), (5, 7), (4, 5), (7, 8), (5, 8)]

    sites = siteinds("S=1/2", N, conserve_qns = false);

    for b in dlattice
        println(b)
    end

    hex_sumOp = OpSum()

    for b in dlattice
        hex_sumOp +=  "Sz", b[1]+1, "Sz", b[2]+1
        hex_sumOp +=  "Sx", b[1]+1, "Sx", b[2]+1
        hex_sumOp +=  "Sy", b[1]+1, "Sy", b[2]+1
    end

    hex_MPO = MPO(hex_sumOp, sites);
    @show hex_MPO;
    @show hex_sumOp;
    @show size(hex_sumOp)[1] == size(dlattice)[1]*3

    if ED_ref
        eig_vals, eig_vec = diag(hex_MPO, N)
        @show eig_vals[1]
    end

    initial_state = [isodd(n) ? "Up" : "Dn" for n=1:N]
    init_MPS = randomMPS(sites, initial_state);


    sweeps = Sweeps(10)
    if is_trunc
        setmaxdim!(sweeps, trunc_th)
    else
        setmaxdim!(sweeps, 10,20,50,100,200,500,800,1000)
    end

    setcutoff!(sweeps, 1e-10)
    setnoise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0)

    energy_DMRG, psi_DMRG = dmrg(hex_MPO, init_MPS, sweeps);
    @show energy_DMRG

    if ED_ref
        @show energy_DMRG - eig_vals[1]
        delta = inner(psi_DMRG', hex_MPO, psi_DMRG) - eig_vals[1]
    end

    full_DM = reduced_rho_matrix(psi_DMRG);
    I_matrix,norm_factor, mi_pairs = MI_diag(full_DM);

    plot(heatmap(z=I_matrix, colorscale = "Viridis"))
    bond_dim = maxlinkdim(psi_DMRG);
    
    if save
        path_N= string("../Image_results/MI_Heiseberg_Hex_3x3_chi_$(bond_dim).png");
        savefig(plot(heatmap(z=I_matrix, colorscale = "Viridis")), path_N)
    end

    
    if is_trunc
        println("Truncation threshold: $(trunc_th)")
    end
    println("SYSTEM:\nHeisenberg Hexagonal lattice $(N) sites\n")
    println("Julia lattice indices")
    for b in dlattice
        println("Coupling:",(b[1]+1,b[2]+1))
    end 
    println("Python lattice indices")
    for b in dlattice
        println("Coupling:",(b[1], b[2]))
    end
    if ED_ref
        println("Exact Diagonalization Energy: $(eig_vals[1])")
    end
    println("DMRG energy: $(energy_DMRG)")
    println("Normalization factor: $(norm_factor)")
    println(mi_pairs)
end
