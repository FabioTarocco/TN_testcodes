""" script_binaryMERA.jl
---------------------------------------------------------------------
Script file demonstrating how the auto-generated code from the "TensorTrace"
software can be implemented in a tensor network algorithm (in this case
an optimization of a binary MERA for the ground state of the 1D transverse
field quantum Ising model on a finite lattice). This script calls the
"binaryMERA.py" function file as automatically generated from the example
TensorTrace project "binaryMERA.ttp" included with the distribution.

by Glen Evenbly (www.glenevenbly.com) - last modified 06/2020
"""

using LinearAlgebra
using Arpack
using Printf
include("binaryMERA.jl");
include("ncon.jl");

# set simulation parameters
chi = 6;
chi_p = 4;
n_levels = 3;
n_iterations = 2000;
n_sites = 3*(2^(n_levels + 1));

# define hamiltonian
sX = [0  1; 1  0];
sY = [0  -im; im  0];
sZ = [-1  0; 0  1];
sI2 = Matrix{Float64}(I, 2, 2);
sI4 = Matrix{Float64}(I, 4, 4);
sI8 = Matrix{Float64}(I, 8, 8);
ham_s = (-kron(sX,sX)+0.5*kron(sI2,sZ)+0.5*kron(sI2,sZ));
ham_init = 0.5*kron(sI8,kron(ham_s,sI2))+kron(sI4,kron(ham_s,sI4))+0.5*kron(sI2,kron(ham_s,sI8));
chi_b = 4;

# initialize tensors
uDis = Array{Any,1}(undef,n_levels);
uDis[1] = reshape(Matrix{Float64}(I, chi_b^2, chi_b^2),chi_b,chi_b,chi_b,chi_b);
wIso = Array{Any,1}(undef,n_levels);
wIso[1] = rand(chi_b,chi_b,chi);
hamThree = Array{Any,1}(undef,n_levels+1);
hamThree[1] = reshape((ham_init - maximum(eigvals(ham_init))*Matrix{Float64}(I, chi_b^3, chi_b^3)),chi_b,chi_b,chi_b,chi_b,chi_b,chi_b);
rhoThree = Array{Any,1}(undef,n_levels+1);
rhoThree[1] = rand(chi_b,chi_b,chi_b,chi_b,chi_b,chi_b);
for k =2:n_levels
    uDis[k] = reshape(Matrix{Float64}(I, chi^2, chi_p^2),chi,chi,chi_p,chi_p);
    wIso[k] = rand(chi_p,chi_p,chi);
    hamThree[k] = rand(chi,chi,chi,chi,chi,chi);
    rhoThree[k] = rand(chi,chi,chi,chi,chi,chi);
end
hamThree[n_levels+1] = rand(chi,chi,chi,chi,chi,chi);
rhoThree[n_levels+1] = rand(chi,chi,chi,chi,chi,chi);

# do optimization iterations
for p = 1:n_iterations
    # sweep over all levels
    for z = 1:n_levels
        if p > 5
            # optimise disentanglers
            tensor_list = Any[uDis[z], wIso[z], hamThree[z], rhoThree[z+1]];
            uEnv = binaryMERA(tensor_list; which_network=1, which_env=1) +
              binaryMERA(tensor_list; which_network=1, which_env=2) +
              binaryMERA(tensor_list; which_network=2, which_env=1) +
              binaryMERA(tensor_list; which_network=2, which_env=2);
            uSize = size(uEnv);
            F = svd(reshape(uEnv,uSize[1]*uSize[2],uSize[3]*uSize[4]));
            uDis[z] = -reshape((F.U*F.Vt),uSize);
        end

        # optimise isometries
        tensor_list = Any[uDis[z], wIso[z], hamThree[z], rhoThree[z+1]];
        Wenv = binaryMERA(tensor_list; which_network=1, which_env=3) +
          binaryMERA(tensor_list; which_network=1, which_env=4) +
          binaryMERA(tensor_list; which_network=1, which_env=5) +
          binaryMERA(tensor_list; which_network=2, which_env=3) +
          binaryMERA(tensor_list; which_network=2, which_env=4) +
          binaryMERA(tensor_list; which_network=2, which_env=5);
        wSize = size(Wenv);
        F = svd(reshape(Wenv,wSize[1]*wSize[2],wSize[3]));
        wIso[z] = -reshape((F.U*F.Vt),wSize);

        # lift Hamiltonian
        tensor_list = Any[uDis[z], wIso[z], hamThree[z], rhoThree[z+1]];
        hamThree[z+1] = binaryMERA(tensor_list; which_network=1, which_env=12) +
          binaryMERA(tensor_list; which_network=2, which_env=12);
    end

    # diagonalize Hamiltonian
    ham_top = reshape(hamThree[n_levels+1]+permutedims(hamThree[n_levels+1],[2,3,1,5,6,4])+permutedims(hamThree[n_levels+1],[3,1,2,6,4,5]),chi^3,chi^3);
    dtemp, vtemp = eigs(0.5*(ham_top + ham_top'); nev=1, tol=1e-10, which=:SR, maxiter = 1000);
    vtemp = vtemp / norm(vtemp);
    rhoThree[n_levels+1] = reshape(vtemp * vtemp',chi,chi,chi,chi,chi,chi);

    # lower the density matrix
    for z = n_levels:-1:1
        tensor_list = Any[uDis[z], wIso[z], hamThree[z], rhoThree[z+1]];
        rhoThree[z] = 0.5*(binaryMERA(tensor_list; which_network=1, which_env=6) +
          binaryMERA(tensor_list; which_network=2, which_env=6));
    end

    # compute energy and magnetization
    Energy_per_site = tr(reshape(rhoThree[1],chi_b^3,chi_b^3) * ham_init)/2;
    ExpectX = tr(reshape(rhoThree[1],chi_b^3,chi_b^3) * kron(Matrix{Float64}(I, 32, 32),sX));
    EnExact = (-2/sin(pi/(2*n_sites))) / n_sites;
    EnError = Energy_per_site - EnExact;
    @printf "Iteration: %d of %d, Energy: %f, Energy Error: %e, XMag: %e\n" p n_iterations Energy_per_site EnError ExpectX
end
