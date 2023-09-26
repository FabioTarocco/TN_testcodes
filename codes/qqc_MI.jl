using ITensors
using ITensorChemistry
using PlotlyJS


function reduced_rhos(psi,i,j)

        orthogonalize!(psi, i)
    
        psi_dag = dag(psi)
        
        if i==j
            rho = prime(psi_dag[i], "Site") * psi[i]
            @assert 2 == length(inds(rho))
            return rho
        end
    
        prime!(psi_dag, "Link")
        if i == 1
            rho = prime(psi_dag[i],"Site")
            rho *= psi[i]
        else 
            link_i_1 = commonind(psi[i],psi[i-1])
            rho = prime(psi_dag[i],"Site");
            rho *= prime(psi[i], link_i_1);
        end
        
        for k in i+1:j-1
            rho *= psi_dag[k]
            rho *= psi[k]
        end

        if j == size(psi)[1]
            rho *= prime(psi_dag[j], "Site")
            rho *= psi[j]
        else
            link_j = commonind(psi[j], psi[j+1])
            rho *= prime(psi_dag[j], "Site")
            rho *= prime(psi[j], link_j)
        end
        
        @assert 4 == length(inds(rho))
        return rho
    end

function rdm_matrix(psi)
    N = length(psi)
    density_matrix = Array{Any}(undef, N,N)
    for i in 1:N, j in i:N
        @show (i,j)
        density_matrix[i,j] = reduced_rhos(psi, i, j)
    end
    return density_matrix
end


function MI_matrix(rdm)
    A,B = size(rdm);
    I_AB = zeros(A,B);
    for a in 1:A, b in a+1:B
        rho_AB = copy(rdm[a,b]);
        delta_AA = delta(dag(inds(rho_AB)[1]), dag(inds(rho_AB)[2]));
        delta_BB = delta(dag(inds(rho_AB)[3]), dag(inds(rho_AB)[4]));

        rho_B = delta_AA*rho_AB;
        rho_A = delta_BB*rho_AB;

        S_AB = -tr(reshape(Array(rho_AB,inds(rho_AB)[1], inds(rho_AB)[3], inds(rho_AB)[2], inds(rho_AB)[4]),16,16) 
                * log(reshape(Array(rho_AB,inds(rho_AB)[1], inds(rho_AB)[3], inds(rho_AB)[2], inds(rho_AB)[4]),16,16)));
        S_A = -tr(Array(rho_B, inds(rho_B)[1], inds(rho_B)[2]) * log(Array(rho_B, inds(rho_B)[1], inds(rho_B)[2])));
        S_B = -tr(Array(rho_A, inds(rho_A)[1], inds(rho_A)[2]) * log(Array(rho_A, inds(rho_A)[1], inds(rho_A)[2])));
        
        I_AB[a,b] = I_AB[b,a] = S_A + S_B - S_AB;
    end
    return I_AB
end

function settings_TN(molecule, basis, type)

    #=
    type = {"Electron", "Qubit", "Fermion"}
    molecule = {"H₂",
                "N₂",
                "H₂O",
                "BeH₂",
                "ammonia",
                "formaldehyde",
                "methane",
                "ethanol",
                "benzene"}
    =#

    molecule = Molecule(molecule)
    basis = basis
    
    @show molecule
    @show basis
    
    println("\nRunning Hartree-Fock")
    hf = @time molecular_orbital_hamiltonian(molecule; basis, sitetype=type)
    hamiltonian = hf.hamiltonian;
    hf_state = hf.hartree_fock_state;
    hf_energy = hf.hartree_fock_energy;
    terms = length(hamiltonian);
    println("N terms: $terms")
    println("Hartree-Fock complete")
    
    println("Basis set size = ", length(hf_state))
    s = siteinds(type, length(hf_state), conserve_qns = true);
    println("\nConstruct MPO")
    H = @time MPO(hamiltonian, s);
    println("MPO Constructed")
    
    @show maxlinkdim(H)
    psi_hf = MPS(s, hf_state);
    return hamiltonian, H, psi_hf, hf_energy
end



ham , H_mpo, psi_hf, hf_energy = settings_TN("H₂O", "sto-3g", "Electron");

nsweeps = 10;
maxdim = [20,60,100,200,400,800];
cutoff = [1E-10];

dmrg_en, dmrg_psi = dmrg(H_mpo, psi_hf; nsweeps, maxdim, cutoff)
@show maxlinkdim(dmrg_psi)

full_DM = rdm_matrix(dmrg_psi);


I_matrix = MI_matrix(full_DM);
