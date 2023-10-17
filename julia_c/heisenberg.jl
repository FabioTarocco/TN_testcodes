#=
HEISEBERG 2D (3x3 lattice, S=1/2) N = 9 sites
    -Compute MPO
    -Initialize MPS
    -DMRG with X = 8 multiple sweeps, X = 8 single sweep and X = [10,50,100,400,800]
    -Compute all 2-sites reduced density matrix
    -Compute MI for each couple 
=#

# Define the adjacency matrix A
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

# Get the number of nodes
num_nodes = size(A, 1)

# Initialize a list to store the pairs of edges
edge_pairs = []

# Iterate over each node pair (i, j) to find edges
for i in 1:num_nodes
    for j in 1:num_nodes
        if A[i, j] == 1
            push!(edge_pairs, (i, j))
        end
    end
end

# Print the list of edge pairs
for pair in edge_pairs
    println(pair)
end

using ITensors
using ITensorUnicodePlots


#---------Reduced density matrix -------------------

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
 
    #Move the orthogonality center to the i-th site
    orthogonalize!(psi, i)

    #create the psi^dag
    psi_dag = dag(psi)
    
    if i==j
        rho = prime(psi_dag[i], "Site") * psi[i]
        @assert 2 == length(inds(rho))
        return rho
    end

    #prime all the virtual index of psi^dag
    prime!(psi_dag, "Link")

    #if i is the first site, then we just need to compute the product state of the first site of psi^dag and psi
    if i == 1
        #prime the physical index of psi[i]^dag, otherwise it will be contracted with the physical index of psi[i]
        rho = prime(psi_dag[i],"Site")
        #add the site i to rho
        rho *= psi[i]
    else #if i is not 1, then we need to get link_i-1
        link_i_1 = commonind(psi[i],psi[i-1])
        #add to rho the i-th site of psi^dag and prime the physical index
        rho = prime(psi_dag[i],"Site");
        #add the i-th site of psi by priming link_i-1
        rho *= prime(psi[i], link_i_1);
    end
    #add to rho all the intermediate sites of psi and psi^dag between site i and j
    for k in i+1:j-1
        rho *= psi_dag[k]
        rho *= psi[k]
    end

    #if j is equals to N (last mps site) then we just need to add psi[j]^dag priming the physical index and psi[j]
    if j == size(psi)[1]
        rho *= prime(psi_dag[j], "Site")
        rho *= psi[j]
    else
        #if j is not the last site, we need to get the link index between the j-th site and the j+1-th site
        link_j = commonind(psi[j], psi[j+1])
        #add to rho the j-th site of psi^dag priming the physical index and then add the j-th site of psi
        rho *= prime(psi_dag[j], "Site")
        rho *= prime(psi[j], link_j)
    end

    @assert 4 == length(inds(rho))
    return rho
end




#-------- Full ALL-PAIRS RDM -----------
function reduced_rho_matrix(psi)
    density_matrix = Array{Any}(undef, N,N)
    
    for i in 1:N
        for j in i:N
            density_matrix[i,j] = reduced_rho(psi, i, j);
        end
    end
    return density_matrix
end


#----------- Mutual info----------------
function MI(rdm)
    
    I_AB = zeros(N,N)
    A,B = size(rdm)
    for a in 1:A
        for b in a+1:B
            rho_AB = copy(rdm[a,b]);
            @show (a,b);
            delta_AA = delta(dag(inds(rho_AB)[2]), dag(inds(rho_AB)[1]));
            delta_BB = delta(dag(inds(rho_AB)[4]), dag(inds(rho_AB)[3]));
            rho_B = delta_AA*rho_AB;
            rho_A = delta_BB*rho_AB;
            S_AB = -tr(reshape(Array(rho_AB,inds(rho_AB)[1], inds(rho_AB)[3], inds(rho_AB)[2], inds(rho_AB)[4]),4,4) * log(reshape(Array(rho_AB,inds(rho_AB)[1], inds(rho_AB)[3], inds(rho_AB)[2], inds(rho_AB)[4]),4,4)));
            S_A = -tr(Array(rho_B, inds(rho_B)[1], inds(rho_B)[2]) * log(Array(rho_B, inds(rho_B)[1], inds(rho_B)[2])));
            S_B = -tr(Array(rho_A, inds(rho_A)[1], inds(rho_A)[2]) * log(Array(rho_A, inds(rho_A)[1], inds(rho_A)[2])));
            @show S_A, S_B, S_AB, (a,b)
            I_AB[a,b] = I_AB[b,a]= S_A + S_B - S_AB;
            
        end
    end 
    return I_AB;
end

function MI_diag(rdm)
    
    I_AB = zeros(N,N)
    A,B = size(rdm)
    for a in 1:A
        for b in a+1:B
            rho_AB = copy(rdm[a,b]);

            L_inds = (dag(inds(rho_AB)[1]), dag(inds(rho_AB)[3]));
            R_inds = (dag(inds(rho_AB)[2]), dag(inds(rho_AB)[4]));
            delta_AA = delta(dag(inds(rho_AB)[2]), dag(inds(rho_AB)[1]));
            delta_BB = delta(dag(inds(rho_AB)[4]), dag(inds(rho_AB)[3]));

            rho_B = delta_AA*rho_AB;

            rho_A = delta_BB*rho_AB;


            D_AB, _ = eigen(rho_AB, L_inds, R_inds);

            S_AB = 0.0;
            for n=1:dim(D_AB, 1)
                p = D_AB[n,n]
                S_AB -= p * log(p)
            end

            D_A, _ = eigen(rho_A)
            S_A = 0.0
            for n=1:dim(D_A, 1)
                p = D_A[n,n]
                S_A -= p * log(p)
            end

            D_B, _ = eigen(rho_B)
            S_B = 0.0
            for n=1:dim(D_B, 1)
                p = D_B[n,n]
                S_B -= p * log(p)
            end

            println("--------------------------\nRDM\n----------------------")
            @show rho_AB;
            @show rho_A;
            @show rho_B;
            println("---------------------\nCheck norm\n------------------------")
            @show rho_AB*delta_AA*delta_BB;
            println("-----------------------\nEigval\n-----------------------")
            @show D_AB;
            @show D_A;
            @show D_B;

            #S_AB = -tr(reshape(Array(rho_AB,inds(rho_AB)[1], inds(rho_AB)[3], inds(rho_AB)[2], inds(rho_AB)[4]),4,4) * log(reshape(Array(rho_AB,inds(rho_AB)[1], inds(rho_AB)[3], inds(rho_AB)[2], inds(rho_AB)[4]),4,4)))
            #S_A = -tr(Array(rho_B, inds(rho_B)[1], inds(rho_B)[2]) * log(Array(rho_B, inds(rho_B)[1], inds(rho_B)[2])))
            #S_B = -tr(Array(rho_A, inds(rho_A)[1], inds(rho_A)[2]) * log(Array(rho_A, inds(rho_A)[1], inds(rho_A)[2])))
            @show S_A, S_B, S_AB, (a,b);

            I_AB[a,b] =I_AB[b,a]= S_A + S_B - S_AB
            
        end
    end 
    return I_AB
end



#----------------MAIN------------------
using ITensors
Nx =3;
Ny = 3;
N = Nx * Ny;

mps_mode = 0 #0=std mps, 1=snake
lattice_type = 0 # 0=open boundary, 1=periodic on Y
heisenberg_mode = 0 #0=ZZXXYY, 1=+--+ZZ
diag = false

sites = siteinds("S=1/2", N, conserve_qns = false)
#sites = siteind("Qubit", N, conserve_qns = false)

if lattice_type == 0
    dlattice = square_lattice(Nx, Ny; yperiodic=false)
elseif lattice_type == 1
    dlattice = square_lattice(Nx, Ny; yperiodic=true)
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

if diag == true
    using LinearAlgebra
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

initial_state = [isodd(n) ? "Up" : "Dn" for n=1:N]
psi_init_8 = randomMPS(sites, initial_state);


sweeps = Sweeps(10)
setmaxdim!(sweeps, 8,16,50,100,200,400,1000,2000)
setcutoff!(sweeps, 1e-10)
setnoise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0)
energy_8, psi_8 = dmrg(heisenberg_2D_H, psi_init_8, sweeps);


if maxlinkdim(psi_8) > 8
    psi_8_trunc = psi_8
    truncate!(psi_8_trunc, maxlinkdim=8)
end

@show energy_8 - eig_vals[1]
@show inner(psi_8_trunc, heisenberg_2D_H, psi_8_trunc) - eig_vals[1]


full_DM_trunc = reduced_rho_matrix(psi_8_trunc);
full_DM = reduced_rho_matrix(psi_8)

I_matrix_trunc = MI(full_DM_trunc);
I_matrix = MI(full_DM);
I_diag_matrix = MI_diag(full_DM);
I_diag_matrix = I_diag_matrix ./ findmax(I_diag_matrix)[1];
plot(heatmap(z=I_diag_matrix, colorscale = "Viridis"))

@show I_matrix

using PlotlyJS

maxC_trunc = findmax(I_matrix_trunc)
maxC = findmax(I_matrix)
I_matrixN = I_matrix ./ maxC[1]

I_matrixN_trunc = I_matrix_trunc ./ maxC_trunc[1]
plot(heatmap(z=I_matrixN, colorscale = "Viridis"))

savefig(plot(heatmap(z=I_matrixN, colorscale = "Viridis")), "heiseberg5x5_chi_1567.png")
plot(heatmap(z=I_matrixN_trunc, colorscale = "Viridis"))
savefig(plot(heatmap(z=I_matrixN_trunc, colorscale = "Viridis")), "heiseberg4x4_psi_8_trunc.png")


using ITensors.HDF5
f = h5open("MPS_GS_8chi_Heisenberg3x3.h5","w")
write(f,"Psi-DMRG_chi8",psi_8)
close(f)

using Combinatorics
V = ITensor(1.)
for j=1:N
    V *= (psi_8[j]*state(sites[j],elem[j]))
  end
v = scalar(V) # Up is 0, Dn is 1
psi_8_totalrank = prod(psi_8)
@show psi_8_totalrank
