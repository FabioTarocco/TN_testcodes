#=
HEISEBERG 2D (3x3 lattice, S=1/2) N = 9 sites
    -Compute MPO
    -Initialize MPS
    -DMRG with X = 8 multiple sweeps, X = 8 single sweep and X = [10,50,100,400,800]
    -Compute all 2-sites reduced density matrix
    -Compute MI for each couple 
=#

using ITensors
using ITensorUnicodePlots

function reduced_rho(psi,i,j)
    #=
    #Move the orthogonality center on the i-th site
    orthogonalize!(psi, i)

    #get the dual mps
    psi_dag = dag(psi)

    #prime all the links in the dag
    prime(psi_dag, "Link")

    #get the link connecting site i-1 and site i
    if i ==1
        rho = psi[i]*prime(psi_dag[i], "Site")
    
    else
        link_i_1 = commonind(psi[i],psi[i-1])
        rho = prime(psi[i],link_i_1) * prime(psi_dag[i], "Site")
        @visualize rho
    end

    #all the sites before i-th site are skipped because contract to identity
    #so we can contract dag(psi_i)*psi_i on the link previous site-link
    #during the contraction, it's needed to prime the site i or the index "site" will be contracted between psi and psi_dag
    

    #contract all the tensor between the i-th site and the j-th site
    for k in i+1:j-1
        rho *= psi[k]
        rho *= psi_dag[k]
    end

    if j == size(psi)[1]
        #if its the last site, the j-th is contracted to rho by link_j-1 leaving Sj
        rho *= psi[j]
        #then the dag site is contracted on link_j-1' and leaving Sj'
        rho *= prime(psi_dag[j], "Site")
    else
        
        #get the link connecting site j and site j+1
        link_j = commonind(psi[j],psi[j+1])
        #rho now has Si, Si', link_j-1, link_j-1' as open indices
        #first we contract the j-th site by link_j-1 leaving open Sj and link_j'
        rho *= prime(psi[j], link_j)
        #the we add the j-th site of the psi_dag by contracting index link_j-1' and link_j' and leaving Sj'
        rho *= prime(psi_dag[j], "Site")
    end
    
    return rho
=#
    orthogonalize!(psi, i)
    psi_dag = dag(psi)
    prime!(psi_dag, "Link")


    if i == 1
        nothing
    else
        link_i_1 = commonind(psi[i],psi[i-1])
        rho = prime(psi_dag[i],"Site");
        @show rho
    end
end


test_rdm = reduced_rho(psi_8, 2,4);
Nx = 3;
Ny = 3;
N = Nx * Ny;

sites = siteinds("S=1/2", N; conserve_qns = true);
lattice = square_lattice(Nx, Ny; yperiodic=false);

heisenberg_2D_sumOp = OpSum();

for b in lattice
    heisenberg_2D_sumOp +=       "Sz", b.s1, "Sz", b.s2
    heisenberg_2D_sumOp += 1/2,  "S+", b.s1, "S-", b.s2
    heisenberg_2D_sumOp += 1/2,  "S-", b.s1, "S+", b.s2
end

#MPO of the 2D Heisenberg Hamiltonian
heisenberg_2D_H = MPO(heisenberg_2D_sumOp, sites);

#Initial configuration of the sites of the system
initial_state = [isodd(n) ? "Up" : "Dn" for n=1:N]

psi_init_8 = randomMPS(sites, initial_state, 8)
psi_init_20 = randomMPS(sites, initial_state, 20)

nsweeps = 10
maxdim_single = [8]
maxdim_full = [20,60,100,200,400,800]
cutoff = [1E-10]

energy_8, psi_8 = dmrg(heisenberg_2D_H, psi_init_8; nsweeps, maxdim_single, cutoff)
energy_dmrg, psi_dmrg = dmrg(heisenberg_2D_H, psi_init_20; nsweeps, maxdim_full, cutoff)

@show energy_8
@show energy_dmrg









density_matrix = Array{Any}(undef, N,N)

i = 10
j = 20

for k in i+1:j-1
    @show k
end


a_ = dag(psi_8[1])
a = psi_8[1]
link_i = commonind(psi_8[1], psi_8[2])
prime!(a_, link_i_1)

@show a_*a

commonind(psi_8[1],psi_8[2])==commonind(psi_8[2], psi_8[1])