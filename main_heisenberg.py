
import numpy as np
from utils.input import Input
from tools.ansatz import Ansatz
from tools.psiopsi_ver2 import Psiopsi_copy
from scipy.optimize import minimize
from vqe import VQE
from qiskit_nature.second_q.hamiltonians.lattices import SquareLattice, BoundaryCondition
from qiskit_nature.second_q.hamiltonians.heisenberg_model import HeisenbergModel
from qiskit.quantum_info.operators import SparsePauliOp
from qiskit.circuit.quantumcircuit import Parameter
from qiskit.quantum_info import Statevector
import qiskit.quantum_info as qi
import itertools
import time
import pandas as pd
import sys

import pickle
import csv

import numpy as np



#----------------------------------------------
#-----------DIAGONALIZATION SNIPPET------------
#----------------------------------------------

def diagonalize (op):
    eig_va, eig_vec = np.linalg.eig(sum)
    sort_perm = eig_va.argsort()
    eig_va.sort()
    eig_vec = eig_vec[:, sort_perm]
    Exact_eng = eig_va[0]
    #print("DIAGONALIZE GS")
    #print(eig_vec[:,0])
    eigenstate = eig_vec[:,0]
    wfdict={}
    for ind in range(len(eigenstate)):
        if abs(eigenstate[ind])>0.00005:    # We ignore components of the wavefunction smaller than cutoff
            wfdict[ind]=eigenstate[ind]
    sortedlist=sorted(wfdict.items(), key=lambda x: abs(x[1]), reverse=True)
    if sortedlist[0][1].real<0:  # Always print the wavefunction with HF >0
        sortedlist=[(sortedlist[i][0],-sortedlist[i][1]) for i in range(len(sortedlist))]
    for j in range(len(sortedlist)):
        #print('|{0}>   {1:3.7f}'.format(np.binary_repr(sortedlist[j][0],9),sortedlist[j][1].real))
        pass
    #print('diagonalized energy',Exact_eng)
    #print('eigenvalues',eig_va)
    eig_transf = np.real(eig_va)
    #print("EIGEN VALUE GS:")
    #print(eig_va[0])
    return eig_vec, eig_va




#------------------------------------------------------------------
# SETTING FIXED PARAMETERS + ADDITIONAL ENTANGLER MAP -------------
#------------------------(TEST)------------------------------------
def additional_layers_double(var_form, fixed_parameters, additional_entangler_map):
    num_params = var_form.num_parameters
    circuit_params_list=[]
    #funziona solo con two local
    #circuit_params_list=var_form.ordered_parameters
    circuit_params_list = var_form.parameters
    dict_param = {}        
    for i in range (0,num_params,1):
        dict_param.update({circuit_params_list[i]:fixed_parameters[i]}) 

    var_form.assign_parameters(dict_param, True)  
    var_form.barrier()      


    list_ent = []
    for index in range(0,len(additional_entangler_map)):
        list_ent.append(additional_entangler_map[index])
    # string = str(num_params-1)
    string = '-1'   
    for i in range(0,len(list_ent)):
    
        string = str(int(string) + 1)
        var_form.ry(Parameter(str(string)),0)
        for c,t in list_ent[i]:
            string = str(int(string) + 1)
            var_form.ry(Parameter(str(string)),t)
            var_form.cx(c,t)
        
        string = str(int(string) + 1)
        var_form.barrier()
        var_form.ry(Parameter(str(string)),t)
        list_ent_rev = list_ent[i][::-1]

        for c,t in list_ent_rev:
            string = str(int(string) + 1)
            var_form.ry(Parameter(str(string)),c)
            var_form.cx(c,t)
        var_form.barrier()

    
        string = str(int(string) + 1)
        var_form.ry(Parameter(str(string)),0)
        for c,t in list_ent[i]:
            string = str(int(string) + 1)
            var_form.ry(Parameter(str(string)),t)
            var_form.cx(c,t)
        
        string = str(int(string) + 1)
        var_form.barrier()
        var_form.ry(Parameter(str(string)),t)
        list_ent_rev = list_ent[i][::-1]

        for c,t in list_ent_rev:
            string = str(int(string) + 1)
            var_form.ry(Parameter(str(string)),c)
            var_form.cx(c,t)
            string = str(int(string) + 1)
            var_form.ry(Parameter(str(string)),t)
        
        string = str(int(string) + 1)
        var_form.ry(Parameter(str(string)),0)
        var_form.barrier()
    return()



#------------------------------------------------------------
#---------------------CUSTOM LAYER DOUBLE V-------------------
#-----------------------(hardcoded W shape)-------------------
#------------------------------------------------------------


def add_layer_W(var_form,  additional_entangler_map, fixed_parameters=None):
    if len(fixed_parameters)==0:
        string = str((len(var_form.parameters))-1) 
        print("CIAOOO")
    if len(fixed_parameters)>0:
        
        string = str((len(var_form.parameters))-1) 
        #string = '-1'
        print(len(var_form.parameters))
        print(len(fixed_parameters))
        print(fixed_parameters)
        num_params = var_form.num_parameters
        circuit_params_list=[]

        circuit_params_list = var_form.parameters
        dict_param = {}        
        for i in range (0,num_params,1):
            dict_param.update({circuit_params_list[i]:fixed_parameters[i]}) 
            print("Iter:{}, Parameter List:{}, Value: {}".format(i,circuit_params_list[i], fixed_parameters[i]))

        var_form.assign_parameters(dict_param, True)  
        var_form.barrier()
          


    list_ent = []
    for index in range(0,len(additional_entangler_map)):
        list_ent.append(additional_entangler_map[index])
    #string = '-1'

    for i in range(0,len(list_ent)):
        flatten_ = set(np.array(list_ent[i]).flatten())
        #down ladder first V
        for s in flatten_:
            string = str(int(string)+1)
            var_form.ry(Parameter(str(string)),s)

        for c,t in list_ent[i]:
            var_form.cx(c,t)
        #Up ladder first V
        var_form.barrier()
        for s in flatten_:
            string = str(int(string)+1)
            var_form.ry(Parameter(str(string)),s)
        var_form.barrier()

        list_ent_rev = list_ent[i][::-1]
        for c,t in list_ent_rev:
            var_form.cx(c,t)
        #down ladder second V
        var_form.barrier()
        for s in flatten_:
            string = str(int(string)+1)
            var_form.ry(Parameter(str(string)),s)

        var_form.barrier()
        for c,t in list_ent[i]:
            var_form.cx(c,t)
        #Up ladder second V
        var_form.barrier()
        for s in flatten_:
            string = str(int(string)+1)
            var_form.ry(Parameter(str(string)),s)
        var_form.barrier()

        list_ent_rev = list_ent[i][::-1]
        for c,t in list_ent_rev:
            var_form.cx(c,t)

        var_form.barrier()
        for s in flatten_:
            string = str(int(string)+1)
            var_form.ry(Parameter(str(string)),s)
        var_form.barrier()
        

    return() 



#------------------------------------------
#-----------------LEONARDO ADDITIONAL LAYER--
#--------------------------------------------

def additional_layers(var_form, fixed_parameters, additional_entangler_map):
    num_params = var_form.num_parameters
    circuit_params_list=[]
    
    circuit_params_list = var_form.parameters
    dict_param = {}        
    for i in range (0,num_params,1):
        dict_param.update({circuit_params_list[i]:fixed_parameters[i]}) 

    var_form.assign_parameters(dict_param, True)  
    print(var_form.draw())
    var_form.barrier()      


    list_ent = []
    for index in range(0,len(additional_entangler_map)):
        list_ent.append(additional_entangler_map[index])
    # string = str(num_params-1)
    string = '-1'       
    for i in range(0,len(list_ent)):
        for c,t in list_ent[i]:
            string = str(int(string) + 1)
            var_form.cx(c,t)
            var_form.ry(Parameter(str(string)),c)
        
        string = str(int(string) + 1)
        var_form.ry(Parameter(str(string)),t)
        var_form.barrier()
        list_ent_rev = list_ent[i][::-1]
        for c,t in list_ent_rev:
            string = str(int(string) + 1)
            var_form.cx(c,t)
            var_form.ry(Parameter(str(string)),t)
        string = str(int(string) + 1)
        var_form.ry(Parameter(str(string)),c)
        var_form.barrier()
    return()   


#--------------------------------------------------
#--------MUTUAL INFORMATION SNIPPET----------------
#--------------------------------------------------

def MI_general(eigen, filename_MI, block = 1, num_qubits = None):
    import matplotlib.pyplot as plt
    eigen = 1/np.sqrt(np.dot(eigen,np.conjugate(eigen)))*eigen   #sistema problema di normalizzazione
    if not bool(num_qubits):
        num_qubits = int(np.log2(len(eigen)))
    if type(block) is int:
        block_dim = []
        qubit_block = int(num_qubits/block)
        #print('\n\n Qubits per block ', str(num_qubits/block), '\n\n', flush=True)
        for i in range(qubit_block):
            t = []
            for j in range(block):
                t.append(i*block + j)
            block_dim.append(t)
    elif type(block) is list:
        block_dim = block.copy()
    I = np.zeros((len(block_dim),len(block_dim)),dtype=complex)
    bipartite_entropy = np.zeros((len(block_dim),len(block_dim)),dtype=complex)
    for i in range(len(block_dim)):
        for j in range(i, len(block_dim)):
            qubits_list = list(range(num_qubits))
            to_keep = [block_dim[i], block_dim[j]]
            to_keep = list(set(list(itertools.chain.from_iterable(to_keep)))) #set is for duplicates in case i=j, from_iterable collapses value to single list
            for q in to_keep:
                qubits_list.remove(q)
            if bool(qubits_list):
                rho_bipartite = qi.partial_trace(eigen,qubits_list)
            else:
                rho_bipartite = eigen
            bipartite_entropy[i][j] = qi.entropy(rho_bipartite)#/len(to_keep) #NON SONO SICURO CI VADA!!!!!!!!
    for i in range(len(block_dim)):
        for j in range(i+1, len(block_dim)):
            I[i][j] = (bipartite_entropy[i][i] + bipartite_entropy[j][j] - bipartite_entropy[i][j])#/(len(block_dim[i])+len(block_dim[j]))
            I[j][i] = I[i][j]
    Imax = I.max()        
    somma = np.sum(I.real)
    I = I/Imax
    #title = 'Mutual info for block {}, Imax{}'.format(str(block_dim), str(Imax))
    title = 'normalization:  '+str(Imax.real)
    plt.figure()
    plt.imshow(I.real, vmin=0, vmax=1)
    #plt.show(block=False)
    #plt.show()
    plt.colorbar()
    plt.title(title)
    plt.savefig(str(filename_MI)+'.png')
    #plt.colorbar().remove()
    #print('mutual information matrix',list(I))
    #print('max value of mutual information between two qubits',Imax.real)
    #print('sum of elements of mutual information matrix',somma)      
    return(list(I))




#------------------------------------------
#-------------------MAIN-------------------
#------------------------------------------
        
print("Input initialized-----------")
input =Input('input.nml')


print("Ansatz created-------------")
ansatz = Ansatz(input= input, n_qubits=9, 
                reduce_parameters=True, 
                starting_occupations=0)
#print(ansatz.var_form.draw())
print("Creation lattice------------")
Nx = 3
Ny = 3
lattice = SquareLattice(Nx,Ny, boundary_condition=BoundaryCondition.OPEN)

print("Created Heisenberg Hamiltonian----------")
heisenberg = HeisenbergModel(lattice=lattice, coupling_constants=(0.25, 0.25, 0.25))
heisengber_secondq_op = heisenberg.second_q_op()


#---------------------------------------
#SETTINGS-------------------------------
#----------------------------------------



diag = False
MI = False
vqe_calc = False
extended_vqe = True
full_vqe = False

test = False

job = sys.argv[2]

vqe_numbers = 1 
filename_MI = "MI_diag_heisenberg3x3"
file_path = "../../res_vqe_ladder/vqe{job}_QIDA_{Nx}x{Ny}_{depth}.pkl".format(job=job, Nx=Nx, Ny=Ny, depth = ansatz.depth[0])


#--------------------------------------------
#OPERATOR-----------------------------------------
#--------------------------------------------

print("Creation of single strings FOR ED--------")
ham_dict = dict()
I_N = "I"*(3*3)
sum = 0
for spin_op in heisengber_secondq_op.keys():
    op = spin_op[0]
    temp = list(I_N)
    temp[int(spin_op[2])] = op
    temp[int(spin_op[6])] = op
    temp = "".join(temp)
    ham_dict[temp] = heisengber_secondq_op[spin_op]
    sum = sum + SparsePauliOp.from_list([(temp, 0.25)]).to_matrix()



#--------------------DIAG----------------------
if diag:
    eig_vec, eig_val = diagonalize(sum)
    #print (eig_val[0])
    #print("EIGEN GS:")
    #print(eig_vec[:,0])

#---------------------MI---------------------
if MI:
    I = MI_general(eig_vec[:,0], filename_MI=filename_MI, num_qubits=9)



#------------------------------------------
#---------VQE COMPUTATION------------------
#------------------------------------------
if vqe_calc:

    try:
        with open(file_path, 'rb') as file:
            data = pickle.load(file)
    except FileNotFoundError:
        data = []

    with open(file_path, 'wb') as file:
        pickle.dump(data, file)


    #print("Starting {i} vqe's with depth {depth}:\n".format(i = job, depth = ansatz.depth))
    for i in  range(vqe_numbers):
        #print("Starting %i iteration"%i)
        params = np.random.uniform(size=ansatz.var_form.num_parameters, low=-np.pi, high=np.pi)
        psiopsi = Psiopsi_copy(input=input, operator=ham_dict, ansatz=ansatz)
        vqe = VQE(initial_params=params, psiopsi=psiopsi, vqe_threshold=input.vqe_threshold, optimizer=input.vqe_optimizer, optimizer_shots=input.vqe_max_opt_iters)
        initial_time = time.time()
        vqe.run_vqe()
        end_time = time.time()
        new_entry = {"Job": job, "Circuit":{"Depth":ansatz.depth, "Entangler_map": ansatz.entangler_map, "Rot_Blocks": ansatz.rotation_blocks}, "Energy": vqe.converged_energy, "Optimal_params": vqe.optimal_parameters[-1], "Time_req": end_time - initial_time}
        data.append(new_entry)
        #print("Completed job: %i "%job)
        #print("-Converged energy %f"%vqe.converged_energy)
        #print("Ended %i iteration\n"%job)

        with open(file_path, 'wb') as file:
            pickle.dump(data, file)



#------------------------------------------
#--------------EXTENDED VQE-----------------
#------------------------------------------





if extended_vqe:
    #-----------------------
    #PATH di dove andare a prendere i parametri ottimali
    #---------------------------
    vqe_df = pd.read_pickle("../../res_vqe_ladder/QIDA_res/3l/vqe11_QIDA_3x3_3.pkl".format(job))[0]
    
    #vqe_df = pd.read_pickle("../../res_vqe_ladder/QIDA_res/3l_2ndL/vqe11_QIDA_V11_3x3_3.pkl".format(job))[0]
    best_params = vqe_df["Optimal_params"]
    #initial_params = vqe_df["Initial_params"]

    circuit=ansatz.var_form.assign_parameters(best_params, inplace=False)
    state=Statevector(circuit)
    s=np.array(state)
    sT=s.conjugate()
    #print(best_params)
    print(vqe_df["Energy"])
    print("ENERGY OBTAINED WITH <psi|O|psi>: ", np.real(sT@sum@s))
    additional_map = [[[1,4],[3,4],[4,5],[4,7]]]
    #additional_map = [[[0,1],[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,8]]]
    print(ansatz.var_form.parameters)
    add_layer_W(ansatz.var_form, fixed_parameters=best_params,additional_entangler_map=additional_map)
    print(ansatz.var_form.parameters)

    print("CIRCUIT")
    print(circuit.draw())
    print("VAR_FORM")
    print(ansatz.var_form.draw())
    print(best_params)
        
    print(ansatz.var_form.draw())
    #file_path = "../../res_vqe_ladder/TEST.pkl"
    file_path = "../../res_vqe_ladder/vqe{job}_QIDA_V{job}_2_{Nx}x{Ny}_{depth}.pkl".format(job=job, Nx=Nx, Ny=Ny, depth = ansatz.depth[0])
    try:
        with open(file_path, 'rb') as file:
            data = pickle.load(file)
    except FileNotFoundError:
        data = []

    with open(file_path, 'wb') as file:
        pickle.dump(data, file)


    #print("Starting {i} Extende vqe's with depth {depth} + VSHAPE :\n".format(i = vqe_numbers, depth = ansatz.depth))
    for i in  range(vqe_numbers):
        #print("Starting %i iteration"%i)
        params = np.zeros((ansatz.var_form.num_parameters))
        psiopsi_extended = Psiopsi_copy(input=input, operator=ham_dict, ansatz=ansatz)
        

        vqe_extended = VQE(initial_params=params, psiopsi=psiopsi_extended, vqe_threshold=input.vqe_threshold, optimizer=input.vqe_optimizer, optimizer_shots=input.vqe_max_opt_iters)
        initial_time = time.time()
        vqe_extended.run_vqe()
        end_time = time.time()
        new_entry = {"Job": job, "Circuit": {"Depth":ansatz.depth, "Entangler_map": ansatz.entangler_map, "Add_layer":additional_map, "Rot_Blocks": ansatz.rotation_blocks},
                        "Energy": vqe_extended.converged_energy, "Pre_Layer_Energy": vqe_df["Energy"],
                        "Initial_params":best_params, "Optimal_params": vqe_extended.optimal_parameters[-1],
                        "Time_req": end_time - initial_time}
        data.append(new_entry)
        #print("Completed %i-th iteration"%i)
        #print("-Converged energy %f"%vqe_extended.converged_energy)
        #print("Ended %i iteration\n"%i)

        with open(file_path, 'wb') as file:
            pickle.dump(data, file)





if full_vqe:

    ansatz = Ansatz(input= input, n_qubits=9, 
                reduce_parameters=True, 
                starting_occupations=0)
    vqe_df = pd.read_pickle("../../res_vqe_ladder/QIDA_res/3l_2ndL/vqe11_QIDA_V11_3x3_3.pkl".format(job))[0]
    #best_params = vqe_df["Initial_params"] + vqe_df["Optimal_params"]
    print(len(vqe_df["Initial_params"]), len(vqe_df["Optimal_params"])) 
    print(ansatz.var_form.draw())

    additional_map = [[[1,4],[3,4],[4,5],[4,7]]]
    add_layer_W(ansatz.var_form, fixed_parameters=None, additional_entangler_map=additional_map)
    print(ansatz.var_form.draw())
    total_params = np.concatenate([vqe_df["Initial_params"],vqe_df["Optimal_params"]])
    print(len(total_params))
    print(vqe_df["Initial_params"])
    print(vqe_df["Optimal_params"])
    print(ansatz.var_form.assign_parameters(total_params))
   
    file_path = "../../res_vqe_ladder/vqe{job}_QIDA_full_{Nx}x{Ny}_{depth}.pkl".format(job=job, Nx=Nx, Ny=Ny, depth = ansatz.depth[0])
    

    try:
        with open(file_path, 'rb') as file:
            data = pickle.load(file)
    except FileNotFoundError:
        data = []

    with open(file_path, 'wb') as file:
        pickle.dump(data, file)


    #print("Starting {i} vqe's with depth {depth}:\n".format(i = job, depth = ansatz.depth))
    for i in  range(vqe_numbers):
        #print("Starting %i iteration"%i)
        params = np.random.uniform(size=ansatz.var_form.num_parameters, low=-np.pi, high=np.pi)
        psiopsi = Psiopsi_copy(input=input, operator=ham_dict, ansatz=ansatz)
        vqe = VQE(initial_params=total_params, psiopsi=psiopsi, vqe_threshold=input.vqe_threshold, optimizer=input.vqe_optimizer, optimizer_shots=input.vqe_max_opt_iters)
        initial_time = time.time()
        vqe.run_vqe()
        end_time = time.time()
        new_entry = {"Job": job, "Circuit":{"Depth":ansatz.depth, "Entangler_map": ansatz.entangler_map, "Rot_Blocks": ansatz.rotation_blocks}, "Energy": vqe.converged_energy, "Optimal_params": vqe.optimal_parameters[-1], "Time_req": end_time - initial_time}
        data.append(new_entry)
        #print("Completed job: %i "%job)
        #print("-Converged energy %f"%vqe.converged_energy)
        #print("Ended %i iteration\n"%job)

        with open(file_path, 'wb') as file:
            pickle.dump(data, file)    






def create_parameters_dictionary(var_form, parameters):
    circuit_params_list = var_form.parameters
    dict_param = {}        
    for i in range (0,len(var_form.parameters),1):
        dict_param.update({circuit_params_list[i]:parameters[i]}) 
        print("Iter:{}, Parameter List:{}, Value: {}".format(i,circuit_params_list[i],parameters[i]))
    return dict_param


if test:
    vqe_df = pd.read_pickle("../../res_vqe_ladder/QIDA_res/3l_2ndL/vqe11_QIDA_V11_3x3_3.pkl".format(job))[0]
    optimal_params = vqe_df["Optimal_params"]
    initial_params = vqe_df["Initial_params"]
    print(initial_params)
    print ("LEN OPTIMAL PARAMS: {}".format(len(optimal_params)))
    print ("LEN INITIAL PARAMS: {}".format(len(initial_params)))
    """circuit=ansatz.var_form.assign_parameters(initial_params, inplace=False)
    state=Statevector(circuit)
    s=np.array(state)
    sT=s.conjugate()
    print("Pre_Layer_Energy: {}".format(vqe_df["Pre_Layer_Energy"]))
    print("ENERGY OBTAINED WITH <psi|O|psi>: ", np.real(sT@sum@s))

    additional_map = [[[1,4],[3,4],[4,5],[4,7]]]
    #additional_map = [[[0,1],[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,8]]]
    print("FIRST LAYER: {}".format(ansatz.var_form.parameters))
    add_layer_W(ansatz.var_form, fixed_parameters=initial_params,additional_entangler_map=additional_map)
    print("SECOND LAYER: {}".format(ansatz.var_form.parameters))

    print("CIRCUIT")
    print(circuit.draw())
    print("VAR_FORM")
    print(ansatz.var_form.draw())

    print("\n\n\tOOOOOOO\n\n\n")"""
    print("COMPOSED FULL CIRCUIT")
    ansatz_2 = Ansatz(input= input, n_qubits=9, reduce_parameters=True, starting_occupations=0)
    vqe_df = pd.read_pickle("../../res_vqe_ladder/QIDA_res/3l_2ndL/vqe11_QIDA_V11_3x3_3.pkl".format(job))[0]
    additional_map = [[[1,4],[3,4],[4,5],[4,7]]]
    add_layer_W(ansatz_2.var_form, fixed_parameters= initial_params, additional_entangler_map=additional_map)
    print()
    ansatz_2.var_form.assign_parameters(optimal_params)
    print(ansatz_2.var_form.draw())
    """total_params = np.concatenate([vqe_df["Initial_params"],vqe_df["Optimal_params"]])
    print(ansatz_2.var_form.draw())
    print("TOTAL PARAMS: {},{}".format(total_params, len(total_params)))
    print(vqe_df["Initial_params"])
    print(vqe_df["Optimal_params"])
    total_dict=create_parameters_dictionary(ansatz_2.var_form, total_params)
    ansatz_2=ansatz_2.var_form.assign_parameters(total_params)


    state=Statevector(ansatz_2)
    s=np.array(state)
    sT=s.conjugate()
    print(np.real(sT@sum@s))
    print("Energy: {}".format(vqe_df["Energy"]))
   """


