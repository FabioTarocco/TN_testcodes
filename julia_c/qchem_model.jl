using ITensors
using ITensorChemistry
using PlotlyJS

include("mutual_info.jl")
include("rdm.jl")


function settings_TN(molecule, basis, type)

    #=
    type = {"Electron", "Qubit", "Fermion"}
    molecule = {"H₂",
                "N₂",
                "NH₃",
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



ham , H_mpo, psi_hf, hf_energy = settings_TN("LiH", "sto-3g", "Fermion");
println("-----------------------\nTensor network settings\n-----------------------");
@show ham;
@show H_mpo;
@show psi_hf;
@show size(psi_hf);
@show hf_energy;


sweeps = Sweeps(10)
setmaxdim!(sweeps, 8,16,50,100,200,400,1000,2000)
setcutoff!(sweeps, 1e-10)
setnoise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0)

dmrg_en, dmrg_psi = dmrg(H_mpo, psi_hf, sweeps);
@show maxlinkdim(dmrg_psi);

full_DM = reduced_rho_matrix(dmrg_psi);
I_matrix = MI_diag(full_DM);
maxC = findmax(I_matrix);
I_matrixN = I_matrix ./ maxC[1];

using PlotlyJS
plot(heatmap(z=I_matrixN, colorscale = "Viridis"))
savefig(plot(heatmap(z=I_matrixN, colorscale = "Viridis")), "./LiH_sto-3g.png")
dmrg_en


s =reorder_MI(I_matrixN)

plot(heatmap(z=s, colorscale = "Viridis"))
savefig(plot(heatmap(z=s, colorscale = "Viridis")), "./LiH_sto-3g_ord.png")


d = plot_MI_coupling(I_matrixN,size(dmrg_psi)[1], 0.0)