using ITensors
using ITensorChemistry
using molecular_orbital_hamiltonian
using PlotlyJS

include("rdm.jl")
include("mutual_info.jl")

#--------------MAIN-----------------


mol = pyscf.gto.