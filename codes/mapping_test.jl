
#------------CHECK MAPPINGS----------

using ITensors
using ITensorChemistry

molecule = Molecule("H₂O")
basis = "sto-3g"

#Electron
hf = molecular_orbital_hamiltonian(molecule; basis) #default is electronic
hamiltonian = hf.hamiltonian
@show length(hamiltonian)
hartree_fock_state = hf.hartree_fock_state
hartree_fock_energy = hf.hartree_fock_energy
electron_hilbert = siteinds("Electron", length(hartree_fock_state); conserve_qns=true)
He = MPO(hamiltonian, electron_hilbert)
ψ₀e = MPS(electron_hilbert, hartree_fock_state)
inner(ψ₀e', He, ψ₀e) ≈ hartree_fock_energy
electronE = copy(hartree_fock_energy)

#Fermion
hf = molecular_orbital_hamiltonian(molecule; basis, sitetype="Fermion")
hamiltonian = hf.hamiltonian
@show length(hamiltonian)
hartree_fock_state = hf.hartree_fock_state
hartree_fock_energy = hf.hartree_fock_energy
fermion_hilbert = siteinds("Fermion", length(hartree_fock_state); conserve_qns=true)
Hf = MPO(hamiltonian, fermion_hilbert)
ψ₀f = MPS(fermion_hilbert, hartree_fock_state)
inner(ψ₀f', Hf, ψ₀f) ≈ electronE


#Qubit
hf = molecular_orbital_hamiltonian(molecule; basis, sitetype="Qubit")
hamiltonian = hf.hamiltonian
@show length(hamiltonian)
hartree_fock_state = hf.hartree_fock_state
hartree_fock_energy = hf.hartree_fock_energy
qubit_hilbert = siteinds("Qubit", length(hartree_fock_state); conserve_qns=true)
Hq = MPO(hamiltonian, qubit_hilbert)
ψ₀q = MPS(qubit_hilbert, hartree_fock_state)
inner(ψ₀q', Hq, ψ₀q) ≈ electronE

#Confronto drmg e check bond dimension
nsweeps = 10;
maxdim = [20,60,100,200,400,800];
cutoff = [1E-10];

sweeps = Sweeps(1)

Ee, psi_e = dmrg(He, ψ₀e; nsweeps, maxdim, cutoff, outputlevel=0)
Ef, psi_f = dmrg(Hf, ψ₀f; nsweeps, maxdim, cutoff, outputlevel=0)
Eq, psi_q = dmrg(Hq, ψ₀q; nsweeps, maxdim, cutoff, outputlevel=0)

@show maxlinkdim(psi_e)
@show maxlinkdim(psi_f)
@show maxlinkdim(psi_q)

Ee ≈ Ef
Ee ≈ Eq