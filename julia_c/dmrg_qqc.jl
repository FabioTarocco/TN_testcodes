using ITensors
using ITensorChemistry


molecule_h2o = Molecule("H₂O")
basis = "sto-3g"

@show molecule_h2o
@show basis

println("\nRunning Hartree-Fock")
hf = @time molecular_orbital_hamiltonian(molecule_h2o; basis, sitetype="Electron")
hamiltonian = hf.hamiltonian
hf_state = hf.hartree_fock_state
hf_energy = hf.hartree_fock_energy
terms = length(hamiltonian)
println("N terms: $terms")
println("Hartree-Fock complete")

println("Basis set size = ", length(hf_state))
#=I need to build a vector of site long as the number of basis set functions
In this case, the correspondence is 1 to 1, because each orbital is a spin type
one, so it can be set in 4 possible states {empy, up, dn, updn}
=#
s = siteinds("Electron", length(hf_state))
println("\nConstruct MPO")
H = @time MPO(hamiltonian, s)
println("MPO Constructed")

@show maxlinkdim(H)
psi_hf = MPS(s, hf_state)

@show inner(psi_hf', H, psi_hf)

@show hf_energy
@show (hf_energy - inner(psi_hf', H, psi_hf))

dmrg_params = (nsweeps = 10, maxdim = [100,200], cutoff = 1e-6)
println("\nRunning DMRG")

@show dmrg_params
e_dmrg, psi_dmrg = @time dmrg(H, psi_hf; dmrg_params... ) 
println("DMRG Completed!")
@show (hf_energy - e_dmrg)
@show (inner(psi_hf', H, psi_hf) - inner(psi_dmrg', H, psi_dmrg))


#=
Check on the dimension and sitetypes
=#
@show molecule

hf = molecular_orbital_hamiltonian(molecule; basis)
hamiltonian = hf.hamiltonian
hartree_fock_state = hf.hartree_fock_state
hartree_fock_energy = hf.hartree_fock_energy

println("Number of orbitals = ", length(hartree_fock_state))
println("Number of fermionic operators = ", length(hamiltonian))
println("|HF⟩ = |", prod(string.(hartree_fock_state)), "⟩")

qubit_state = jordanwigner(hartree_fock_state);
qubit_hamiltonian = @time jordanwigner(hamiltonian);
println("Number of qubit operators = ", length(qubit_hamiltonian));
println("|HF⟩ = |", prod(string.(qubit_state)), "⟩");

#=
----------------- TEST DOUBLE TYPE---------------------
Electron Hamiltonian {vac, up, dn, updn} for each sites
=#

molecule = Molecule("H₂")
basis = "sto-3g"

hf = molecular_orbital_hamiltonian(molecule; basis)
#OpSum of the Electronic Hamiltonian
hamiltonian = hf.hamiltonian
hartree_fock_state = hf.hartree_fock_state
hartree_fock_energy = hf.hartree_fock_energy
nterms = length(hamiltonian)
println("Number of terms in the OpSum: $nterms")

s = siteinds("Electron", length(hartree_fock_state); conserve_qns=true)
println("\n$s",s)
H = MPO(hamiltonian, s)
@show maxlinkdim(H)
ψhf = MPS(s, hartree_fock_state)

inner(ψhf', H, ψhf) ≈ hartree_fock_energy
@show inner(ψhf', H, ψhf) - hartree_fock_energy
sweeps = Sweeps(10)
setmaxdim!(sweeps, 100, 200)
setcutoff!(sweeps, 1e-6)
setnoise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0)
energy_dmrg_electronic, psi_dmrg_electronic= dmrg(H, ψhf, sweeps; outputlevel=0)
@show energy_dmrg_electronic < hartree_fock_energy

#=
Fermionic Hamiltonian {vac, occupied} for 2N sites
=#


molecule = Molecule("H₂")
basis = "sto-3g"
hf = molecular_orbital_hamiltonian(molecule; basis)
hamiltonian = hf.hamiltonian
nterms_e = length(hamiltonian)
hartree_fock_state = hf.hartree_fock_state
hartree_fock_energy = hf.hartree_fock_energy

se = siteinds("Electron", length(hartree_fock_state); conserve_qns=false)
He = MPO(hamiltonian, se)
ψe = MPS(se, hartree_fock_state)

inner(ψe', He, ψe) ≈ hartree_fock_energy
@show inner(ψe', He, ψe) - hartree_fock_energy


#HERE I select the sitetype of the system
hf = molecular_orbital_hamiltonian(molecule; basis, sitetype="Fermion")
hamiltonian = hf.hamiltonian
nterms_f = length(hamiltonian)
hartree_fock_state = hf.hartree_fock_state
hartree_fock_energy = hf.hartree_fock_energy

sf = siteinds("Fermion", length(hartree_fock_state); conserve_qns=true)
Hf = MPO(hamiltonian, sf)
@show maxlinkdim(Hf)
ψf = MPS(sf, hartree_fock_state)

inner(ψf', Hf, ψf) ≈ hartree_fock_energy
@show inner(ψf', Hf, ψf) - hartree_fock_energy

sweeps = Sweeps(10)
setmaxdim!(sweeps, 100, 200)
setcutoff!(sweeps, 1e-10)
setnoise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0)
Ee, ψe = dmrg(He, randomMPS(se), sweeps; outputlevel=0)

sweeps = Sweeps(10)
setmaxdim!(sweeps, 100, 200)
setcutoff!(sweeps, 1e-10)
setnoise!(sweeps, 1e-6, 1e-7, 1e-8, 0.0)
Ef, ψf = dmrg(Hf, ψf, sweeps; outputlevel=0)
Ee ≈ Ef

