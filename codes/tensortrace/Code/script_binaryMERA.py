"""
script_binaryMERA.py
---------------------------------------------------------------------
Script file demonstrating how the auto-generated code from the "TensorTrace"
software can be implemented in a tensor network algorithm (in this case
an optimization of a binary MERA for the ground state of the 1D transverse
field quantum Ising model on a finite lattice). This script calls the
"binaryMERA.py" function file as automatically generated from the example
TensorTrace project "binaryMERA.ttp" included with the distribution.

by Glen Evenbly (www.glenevenbly.com) - last modified 06/2020
"""

from scipy.sparse.linalg import eigsh
from numpy import linalg as LA
import numpy as np
from binaryMERA import binaryMERA

# set simulation parameters
chi = 6
chi_p = 4
n_levels = 3
n_iterations = 2000
n_sites = 3 * (2**(n_levels + 1))

# define hamiltonian and do 2-to-1 blocking
chi_b = 4
sX = np.array([[0, 1], [1, 0]])
sZ = np.array([[1, 0], [0, -1]])
ham_s = (-np.kron(sX, sX) + 0.5 * np.kron(np.eye(2), sZ) +
         0.5 * np.kron(np.eye(2), sZ))
ham_init = (0.5 * np.kron(np.eye(8), np.kron(ham_s, np.eye(2))) +
            np.kron(np.eye(4), np.kron(ham_s, np.eye(4))) +
            0.5 * np.kron(np.eye(2), np.kron(ham_s, np.eye(8))))
bias = max(LA.eigvalsh(ham_init))

# initialize tensors
uDis = [0] * n_levels
wIso = [0] * n_levels
ham = [0] * (n_levels + 1)
rho = [0] * (n_levels + 1)
uDis[0] = np.eye(chi_b**2, chi_b**2).reshape(chi_b, chi_b, chi_b, chi_b)
wIso[0] = np.random.rand(chi_b, chi_b, chi)
ham[0] = (ham_init - bias * np.eye(chi_b**3)).reshape([chi_b] * 6)
rho[0] = np.random.rand(chi_b, chi_b, chi_b, chi_b, chi_b, chi_b)
for k in range(n_levels - 1):
  uDis[k + 1] = (np.eye(chi**2, chi_p**2)).reshape(chi, chi, chi_p, chi_p)
  wIso[k + 1] = np.random.rand(chi_p, chi_p, chi)
  ham[k + 1] = np.random.rand(chi, chi, chi, chi, chi, chi)
  rho[k + 1] = np.random.rand(chi, chi, chi, chi, chi, chi)

ham[n_levels] = np.random.rand(chi, chi, chi, chi, chi, chi)
rho[n_levels] = np.random.rand(chi, chi, chi, chi, chi, chi)

# do iterations of single tensor optimizations
for p in range(n_iterations):

  # sweep over all levels
  for z in range(n_levels):
    if p > 10:
      # optimise disentanglers
      tensor_list = [uDis[z], wIso[z], ham[z], rho[z + 1]]
      uEnv = (binaryMERA(tensor_list, which_network=1, which_env=1) +
              binaryMERA(tensor_list, which_network=1, which_env=2) +
              binaryMERA(tensor_list, which_network=2, which_env=1) +
              binaryMERA(tensor_list, which_network=2, which_env=2))
      uSize = uEnv.shape
      uEnv = uEnv.reshape(uSize[0] * uSize[1], uSize[2] * uSize[3])
      utemp, stemp, vtemph = LA.svd(uEnv, full_matrices=False)
      uDis[z] = (-utemp @ vtemph).reshape(uSize)

    # optimise isometries
    tensor_list = [uDis[z], wIso[z], ham[z], rho[z + 1]]
    wEnv = (binaryMERA(tensor_list, which_network=1, which_env=3) +
            binaryMERA(tensor_list, which_network=1, which_env=4) +
            binaryMERA(tensor_list, which_network=1, which_env=5) +
            binaryMERA(tensor_list, which_network=2, which_env=3) +
            binaryMERA(tensor_list, which_network=2, which_env=4) +
            binaryMERA(tensor_list, which_network=2, which_env=5))
    wSize = wEnv.shape
    wEnv = wEnv.reshape(wSize[0] * wSize[1], wSize[2])
    utemp, stemp, vtemph = LA.svd(wEnv, full_matrices=False)
    wIso[z] = (-utemp @ vtemph).reshape(wSize)

    # lift Hamiltonian
    tensor_list = [uDis[z], wIso[z], ham[z], rho[z + 1]]
    ham[z + 1] = (binaryMERA(tensor_list, which_network=1, which_env=12) +
                  binaryMERA(tensor_list, which_network=2, which_env=12))

  # diagonalize top-level Hamiltonian
  ham_top = (ham[n_levels] +
             ham[n_levels].transpose(1, 2, 0, 4, 5, 3) +
             ham[n_levels].transpose(2, 0, 1, 5, 3, 4)
             ).reshape(chi**3, chi**3)
  dtemp, vtemp = eigsh(0.5 * (ham_top + np.conj(ham_top.T)), k=1, which='SA')
  vtemp = vtemp / LA.norm(vtemp)
  top_shape = [chi] * 6
  rho[n_levels] = (vtemp @ np.conj(vtemp.T)).reshape([chi] * 6)

  # lower the density matrix
  for z in reversed(range(n_levels)):
    tensor_list = [uDis[z], wIso[z], ham[z], rho[z + 1]]
    rho[z] = 0.5 * (binaryMERA(tensor_list, which_network=1, which_env=6) +
                    binaryMERA(tensor_list, which_network=2, which_env=6))

  # compute energy and magnetization
  Energy_per_site = np.trace(rho[0].reshape(chi_b**3, chi_b**3) @ ham_init) / 2
  ExpectX = np.trace(rho[0].reshape(chi_b**3, chi_b**3) @
                     np.kron(np.eye(32), sX))
  EnExact = (-2 / np.sin(np.pi / (2 * n_sites))) / n_sites
  EnError = Energy_per_site - EnExact
  print('Iteration: %d of %d, Energy: %f, Energy Error: %e, XMag: %e\n'
        % (p, n_iterations, Energy_per_site, EnError, ExpectX))
