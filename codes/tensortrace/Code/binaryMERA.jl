# binaryMERA.jl 
include("ncon.jl"); 

function binaryMERA(tensors_in; which_network=1, which_env=0) 
# binaryMERA(tensors_in; which_network=1, which_env=0) 
# 
# Auto-generated network contraction function from 'TensorTrace' software, 
# see (www.tensortrace.com) for details, (c) Glen Evenbly, 2019. 
# Requires the network contraction routine 'ncon' to be in the working directory, 
# (included in the TensorTrace install directory, or also found at 'www.tensortrace.com'). 
# 
# Input variables: 
# ------------------------------- 
# 1st input 'tensors_in' is a cell array of the unique tensors: tensors_in = {u, w, ham, rho}; 
# 2nd input 'which_network' dictates which network to evaluate (1-4) 
# 3rd input 'which_env' allows one to specify a tensor environment to evaluate from a closed tensor network. 
# 	 -set 'which_env = 0' to evaluate the scalar from a closed network (i.e. no environment). 
# 	 -set 'which_env = n' to evaluate environment of the nth tensor from the closed network. 
# 
# General project info: 
# ------------------------------- 
# Generated on: 30/6/2020
# Generated from: C:\Users\gevenbly3\Desktop\TensorTrace\SaveData\
# Index dims: chi = 8,chi_p = 6

# --- Info for Network-1 --- 
# -------------------------- 
# network is CLOSED 
# total contraction cost (in scalar multiplications): 8.86*10^7
# contraction order: (((((((T1*T4)*(T2*(T6*T8)))*T7)*T10)*(T3*T9))*(T5*T11))*T12)
# 1st leading order cost: (chi^4)*(chi_p^5)
# 2nd leading order cost: (chi^3)*(chi_p^6)
# tensors_N1 = Any[u,u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w),rho]; 
# connects_N1 = Any[[1,3,10,11],[4,7,12,13],[8,10,17],[11,12,18],[13,14,19],[2,5,6,3,4,7],[1,2,9,23],[5,6,16,15],[8,9,22],[23,16,21],[15,14,20],[22,21,20,17,18,19]]; 
# dims_N1 = Any[[chi,chi,chi_p,chi_p],[chi,chi,chi_p,chi_p],[chi_p,chi_p,chi],[chi_p,chi_p,chi],[chi_p,chi_p,chi],[chi,chi,chi,chi,chi,chi],[chi,chi,chi_p,chi_p],[chi,chi,chi_p,chi_p],[chi_p,chi_p,chi],[chi_p,chi_p,chi],[chi_p,chi_p,chi],[chi,chi,chi,chi,chi,chi]]; 
# con_order_N1 = [11, 5, 6, 4, 7, 3, 12, 1, 2, 14, 8, 16, 23, 10, 9, 13, 15, 18, 21, 17, 22, 19, 20]; 

# --- Info for Network-2 --- 
# -------------------------- 
# network is CLOSED 
# total contraction cost (in scalar multiplications): 8.86*10^7
# contraction order: (((((((T1*(T6*T7))*(T8*T10))*T2)*T4)*(T3*T9))*(T5*T11))*T12)
# 1st leading order cost: (chi^4)*(chi_p^5)
# 2nd leading order cost: (chi^3)*(chi_p^6)
# tensors_N2 = Any[u,u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w),rho]; 
# connects_N2 = Any[[2,3,9,10],[6,23,11,12],[7,9,16],[10,11,17],[12,13,18],[1,4,5,2,3,6],[1,4,8,22],[5,23,15,14],[7,8,21],[22,15,20],[14,13,19],[21,20,19,16,17,18]]; 
# dims_N2 = Any[[chi,chi,chi_p,chi_p],[chi,chi,chi_p,chi_p],[chi_p,chi_p,chi],[chi_p,chi_p,chi],[chi_p,chi_p,chi],[chi,chi,chi,chi,chi,chi],[chi,chi,chi_p,chi_p],[chi,chi,chi_p,chi_p],[chi_p,chi_p,chi],[chi_p,chi_p,chi],[chi_p,chi_p,chi],[chi,chi,chi,chi,chi,chi]]; 
# con_order_N2 = [15, 1, 4, 2, 3, 5, 22, 6, 23, 7, 10, 11, 9, 8, 13, 14, 12, 20, 17, 16, 21, 18, 19]; 

# --- Info for Network-3 --- 
# -------------------------- 
# -- Network is empty. -- 

# --- Info for Network-4 --- 
# -------------------------- 
# -- Network is empty. -- 


# Contraction Code: 
# ------------------------------- 

	if length(tensors_in) == 0 # auto-generate random tensors of default dims 
		chi = 8; 
		chi_p = 6; 
		u = rand(chi,chi,chi_p,chi_p); 
		w = rand(chi_p,chi_p,chi); 
		ham = rand(chi,chi,chi,chi,chi,chi); 
		rho = rand(chi,chi,chi,chi,chi,chi); 
	else
		u = tensors_in[1]; 
		w = tensors_in[2]; 
		ham = tensors_in[3]; 
		rho = tensors_in[4]; 
	end 

	if which_network == 1
		if which_env == 0 
			# TTv1.0.5.0J$;@<EHE<?H?@TE`ET?`?B09<963BH9T9N3B`9l9f3>HQTQ`QHKTK`K@<WHW<]H]@TW`WT]`]B0c<c6iBHcTcNiB`clcfiF6'N'f'6-N-f-')++))+++''+''''++++(((LV,HQ^LO=]'IQ^:HUYqHQ^VJ$; 
			tensors = Any[u,u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w),rho]; 
			connects = Any[[1,3,10,11],[4,7,12,13],[8,10,17],[11,12,18],[13,14,19],[2,5,6,3,4,7],[1,2,9,23],[5,6,16,15],[8,9,22],[23,16,21],[15,14,20],[22,21,20,17,18,19]]; 
			con_order = [11, 5, 6, 4, 7, 3, 12, 1, 2, 14, 8, 16, 23, 10, 9, 13, 15, 18, 21, 17, 22, 19, 20]; 
			return ncon(tensors,connects,con_order=con_order); 
		elseif which_env == 1 
			tensors = Any[u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w),rho]; 
			connects = Any[[4,7,12,13],[8,-3,17],[-4,12,18],[13,14,19],[2,5,6,-2,4,7],[-1,2,9,23],[5,6,16,15],[8,9,22],[23,16,21],[15,14,20],[22,21,20,17,18,19]]; 
			con_order = [5, 6, 4, 7, 14, 8, 19, 20, 17, 22, 21, 9, 23, 13, 2, 16, 15, 12, 18]; 
			return ncon(tensors,connects,con_order=con_order); 
		elseif which_env == 2 
			tensors = Any[u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w),rho]; 
			connects = Any[[1,3,10,11],[8,10,17],[11,-3,18],[-4,14,19],[2,5,6,3,-1,-2],[1,2,9,23],[5,6,16,15],[8,9,22],[23,16,21],[15,14,20],[22,21,20,17,18,19]]; 
			con_order = [11, 5, 6, 14, 8, 19, 20, 17, 22, 21, 9, 23, 1, 10, 18, 3, 15, 16, 2]; 
			return ncon(tensors,connects,con_order=con_order); 
		elseif which_env == 3 
			tensors = Any[u,u,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w),rho]; 
			connects = Any[[1,3,-2,11],[4,7,12,13],[11,12,18],[13,14,19],[2,5,6,3,4,7],[1,2,9,23],[5,6,16,15],[-1,9,22],[23,16,21],[15,14,20],[22,21,20,-3,18,19]]; 
			con_order = [11, 5, 6, 4, 7, 3, 12, 1, 2, 14, 16, 23, 19, 20, 18, 13, 15, 21, 9, 22]; 
			return ncon(tensors,connects,con_order=con_order); 
		elseif which_env == 4 
			tensors = Any[u,u,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w),rho]; 
			connects = Any[[1,3,10,-1],[4,7,-2,13],[8,10,17],[13,14,19],[2,5,6,3,4,7],[1,2,9,23],[5,6,16,15],[8,9,22],[23,16,21],[15,14,20],[22,21,20,17,-3,19]]; 
			con_order = [5, 6, 4, 7, 14, 8, 19, 20, 17, 22, 21, 9, 23, 13, 2, 16, 15, 1, 3, 10]; 
			return ncon(tensors,connects,con_order=con_order); 
		elseif which_env == 5 
			tensors = Any[u,u,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w),rho]; 
			connects = Any[[1,3,10,11],[4,7,12,-1],[8,10,17],[11,12,18],[2,5,6,3,4,7],[1,2,9,23],[5,6,16,15],[8,9,22],[23,16,21],[15,-2,20],[22,21,20,17,18,-3]]; 
			con_order = [11, 5, 6, 4, 7, 3, 12, 1, 2, 8, 16, 23, 10, 9, 18, 21, 17, 22, 15, 20]; 
			return ncon(tensors,connects,con_order=con_order); 
		elseif which_env == 6 
			tensors = Any[u,u,w,w,w,conj(u),conj(u),conj(w),conj(w),conj(w),rho]; 
			connects = Any[[1,-4,10,11],[-5,-6,12,13],[8,10,17],[11,12,18],[13,14,19],[1,-1,9,23],[-2,-3,16,15],[8,9,22],[23,16,21],[15,14,20],[22,21,20,17,18,19]]; 
			con_order = [11, 14, 8, 19, 20, 17, 22, 21, 9, 23, 1, 10, 18, 12, 13, 15, 16]; 
			return ncon(tensors,connects,con_order=con_order); 
		elseif which_env == 7 
			tensors = Any[u,u,w,w,w,ham,conj(u),conj(w),conj(w),conj(w),rho]; 
			connects = Any[[-1,3,10,11],[4,7,12,13],[8,10,17],[11,12,18],[13,14,19],[-2,5,6,3,4,7],[5,6,16,15],[8,-3,22],[-4,16,21],[15,14,20],[22,21,20,17,18,19]]; 
			con_order = [11, 5, 6, 4, 7, 3, 12, 14, 8, 19, 20, 17, 22, 21, 10, 18, 13, 16, 15]; 
			return ncon(tensors,connects,con_order=con_order); 
		elseif which_env == 8 
			tensors = Any[u,u,w,w,w,ham,conj(u),conj(w),conj(w),conj(w),rho]; 
			connects = Any[[1,3,10,11],[4,7,12,13],[8,10,17],[11,12,18],[13,14,19],[2,-1,-2,3,4,7],[1,2,9,23],[8,9,22],[23,-3,21],[-4,14,20],[22,21,20,17,18,19]]; 
			con_order = [11, 14, 8, 19, 20, 17, 22, 21, 9, 23, 1, 10, 18, 12, 13, 3, 2, 4, 7]; 
			return ncon(tensors,connects,con_order=con_order); 
		elseif which_env == 9 
			tensors = Any[u,u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),rho]; 
			connects = Any[[1,3,10,11],[4,7,12,13],[-1,10,17],[11,12,18],[13,14,19],[2,5,6,3,4,7],[1,2,-2,23],[5,6,16,15],[23,16,21],[15,14,20],[-3,21,20,17,18,19]]; 
			con_order = [11, 5, 6, 4, 7, 3, 12, 1, 2, 14, 16, 23, 19, 20, 18, 13, 15, 21, 10, 17]; 
			return ncon(tensors,connects,con_order=con_order); 
		elseif which_env == 10 
			tensors = Any[u,u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),rho]; 
			connects = Any[[1,3,10,11],[4,7,12,13],[8,10,17],[11,12,18],[13,14,19],[2,5,6,3,4,7],[1,2,9,-1],[5,6,-2,15],[8,9,22],[15,14,20],[22,-3,20,17,18,19]]; 
			con_order = [11, 5, 6, 4, 7, 3, 12, 1, 2, 14, 8, 19, 20, 17, 22, 10, 18, 13, 15, 9]; 
			return ncon(tensors,connects,con_order=con_order); 
		elseif which_env == 11 
			tensors = Any[u,u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),rho]; 
			connects = Any[[1,3,10,11],[4,7,12,13],[8,10,17],[11,12,18],[13,-2,19],[2,5,6,3,4,7],[1,2,9,23],[5,6,16,-1],[8,9,22],[23,16,21],[22,21,-3,17,18,19]]; 
			con_order = [11, 5, 6, 4, 7, 3, 12, 1, 2, 8, 16, 23, 10, 9, 18, 21, 17, 22, 13, 19]; 
			return ncon(tensors,connects,con_order=con_order); 
		elseif which_env == 12 
			tensors = Any[u,u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w)]; 
			connects = Any[[1,3,10,11],[4,7,12,13],[8,10,-4],[11,12,-5],[13,14,-6],[2,5,6,3,4,7],[1,2,9,23],[5,6,16,15],[8,9,-1],[23,16,-2],[15,14,-3]]; 
			con_order = [11, 5, 6, 4, 7, 3, 12, 1, 2, 14, 8, 16, 23, 10, 9, 13, 15]; 
			return ncon(tensors,connects,con_order=con_order); 
		else 
			error("error: requested environment ($(which_env)) is out of range for current network; please set 'which_env' in range [0-12]") 
		end 


	elseif which_network == 2
		if which_env == 0 
			# TTv1.0.5.0J$;@<EHE<?H?@TE`ET?`?B09<963BH9T9N3B`9l9f3><QHQTQ<KHKTK@<WHW<]H]@TW`WT]`]B0c<c6iBHcTcNiB`clcfiF6'N'f'6-N-f-))++)'+++''+''''++++(((LV,HQ^LO=]'IQ^:HUYqHQ^VJ$; 
			tensors = Any[u,u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w),rho]; 
			connects = Any[[2,3,9,10],[6,23,11,12],[7,9,16],[10,11,17],[12,13,18],[1,4,5,2,3,6],[1,4,8,22],[5,23,15,14],[7,8,21],[22,15,20],[14,13,19],[21,20,19,16,17,18]]; 
			con_order = [15, 1, 4, 2, 3, 5, 22, 6, 23, 7, 10, 11, 9, 8, 13, 14, 12, 20, 17, 16, 21, 18, 19]; 
			return ncon(tensors,connects,con_order=con_order); 
		elseif which_env == 1 
			tensors = Any[u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w),rho]; 
			connects = Any[[6,23,11,12],[7,-3,16],[-4,11,17],[12,13,18],[1,4,5,-1,-2,6],[1,4,8,22],[5,23,15,14],[7,8,21],[22,15,20],[14,13,19],[21,20,19,16,17,18]]; 
			con_order = [15, 1, 4, 7, 13, 18, 19, 16, 21, 17, 11, 12, 23, 14, 20, 6, 8, 5, 22]; 
			return ncon(tensors,connects,con_order=con_order); 
		elseif which_env == 2 
			tensors = Any[u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w),rho]; 
			connects = Any[[2,3,9,10],[7,9,16],[10,-3,17],[-4,13,18],[1,4,5,2,3,-1],[1,4,8,22],[5,-2,15,14],[7,8,21],[22,15,20],[14,13,19],[21,20,19,16,17,18]]; 
			con_order = [15, 1, 4, 2, 3, 5, 22, 7, 13, 18, 19, 16, 21, 17, 9, 10, 8, 14, 20]; 
			return ncon(tensors,connects,con_order=con_order); 
		elseif which_env == 3 
			tensors = Any[u,u,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w),rho]; 
			connects = Any[[2,3,-2,10],[6,23,11,12],[10,11,17],[12,13,18],[1,4,5,2,3,6],[1,4,8,22],[5,23,15,14],[-1,8,21],[22,15,20],[14,13,19],[21,20,19,-3,17,18]]; 
			con_order = [15, 1, 4, 2, 3, 5, 22, 6, 23, 10, 11, 13, 18, 19, 14, 20, 12, 17, 8, 21]; 
			return ncon(tensors,connects,con_order=con_order); 
		elseif which_env == 4 
			tensors = Any[u,u,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w),rho]; 
			connects = Any[[2,3,9,-1],[6,23,-2,12],[7,9,16],[12,13,18],[1,4,5,2,3,6],[1,4,8,22],[5,23,15,14],[7,8,21],[22,15,20],[14,13,19],[21,20,19,16,-3,18]]; 
			con_order = [15, 1, 4, 2, 3, 5, 22, 6, 23, 7, 13, 18, 19, 16, 21, 9, 8, 14, 20, 12]; 
			return ncon(tensors,connects,con_order=con_order); 
		elseif which_env == 5 
			tensors = Any[u,u,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w),rho]; 
			connects = Any[[2,3,9,10],[6,23,11,-1],[7,9,16],[10,11,17],[1,4,5,2,3,6],[1,4,8,22],[5,23,15,14],[7,8,21],[22,15,20],[14,-2,19],[21,20,19,16,17,-3]]; 
			con_order = [15, 1, 4, 2, 3, 5, 22, 6, 23, 7, 10, 11, 9, 8, 20, 17, 16, 21, 14, 19]; 
			return ncon(tensors,connects,con_order=con_order); 
		elseif which_env == 6 
			tensors = Any[u,u,w,w,w,conj(u),conj(u),conj(w),conj(w),conj(w),rho]; 
			connects = Any[[-4,-5,9,10],[-6,23,11,12],[7,9,16],[10,11,17],[12,13,18],[-1,-2,8,22],[-3,23,15,14],[7,8,21],[22,15,20],[14,13,19],[21,20,19,16,17,18]]; 
			con_order = [15, 7, 13, 18, 19, 16, 21, 17, 11, 12, 23, 14, 20, 9, 10, 8, 22]; 
			return ncon(tensors,connects,con_order=con_order); 
		elseif which_env == 7 
			tensors = Any[u,u,w,w,w,ham,conj(u),conj(w),conj(w),conj(w),rho]; 
			connects = Any[[2,3,9,10],[6,23,11,12],[7,9,16],[10,11,17],[12,13,18],[-1,-2,5,2,3,6],[5,23,15,14],[7,-3,21],[-4,15,20],[14,13,19],[21,20,19,16,17,18]]; 
			con_order = [15, 7, 13, 18, 19, 16, 21, 17, 11, 12, 23, 14, 20, 9, 10, 2, 3, 6, 5]; 
			return ncon(tensors,connects,con_order=con_order); 
		elseif which_env == 8 
			tensors = Any[u,u,w,w,w,ham,conj(u),conj(w),conj(w),conj(w),rho]; 
			connects = Any[[2,3,9,10],[6,-2,11,12],[7,9,16],[10,11,17],[12,13,18],[1,4,-1,2,3,6],[1,4,8,22],[7,8,21],[22,-3,20],[-4,13,19],[21,20,19,16,17,18]]; 
			con_order = [1, 4, 2, 3, 7, 13, 18, 19, 16, 21, 17, 11, 12, 9, 10, 6, 8, 22, 20]; 
			return ncon(tensors,connects,con_order=con_order); 
		elseif which_env == 9 
			tensors = Any[u,u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),rho]; 
			connects = Any[[2,3,9,10],[6,23,11,12],[-1,9,16],[10,11,17],[12,13,18],[1,4,5,2,3,6],[1,4,-2,22],[5,23,15,14],[22,15,20],[14,13,19],[-3,20,19,16,17,18]]; 
			con_order = [15, 1, 4, 2, 3, 5, 22, 6, 23, 10, 11, 13, 18, 19, 14, 20, 12, 17, 9, 16]; 
			return ncon(tensors,connects,con_order=con_order); 
		elseif which_env == 10 
			tensors = Any[u,u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),rho]; 
			connects = Any[[2,3,9,10],[6,23,11,12],[7,9,16],[10,11,17],[12,13,18],[1,4,5,2,3,6],[1,4,8,-1],[5,23,-2,14],[7,8,21],[14,13,19],[21,-3,19,16,17,18]]; 
			con_order = [1, 4, 2, 3, 7, 13, 18, 19, 16, 21, 17, 11, 12, 9, 10, 6, 8, 5, 23, 14]; 
			return ncon(tensors,connects,con_order=con_order); 
		elseif which_env == 11 
			tensors = Any[u,u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),rho]; 
			connects = Any[[2,3,9,10],[6,23,11,12],[7,9,16],[10,11,17],[12,-2,18],[1,4,5,2,3,6],[1,4,8,22],[5,23,15,-1],[7,8,21],[22,15,20],[21,20,-3,16,17,18]]; 
			con_order = [15, 1, 4, 2, 3, 5, 22, 6, 23, 7, 10, 11, 9, 8, 20, 17, 16, 21, 12, 18]; 
			return ncon(tensors,connects,con_order=con_order); 
		elseif which_env == 12 
			tensors = Any[u,u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w)]; 
			connects = Any[[2,3,9,10],[6,23,11,12],[7,9,-4],[10,11,-5],[12,13,-6],[1,4,5,2,3,6],[1,4,8,22],[5,23,15,14],[7,8,-1],[22,15,-2],[14,13,-3]]; 
			con_order = [15, 1, 4, 2, 3, 5, 22, 6, 23, 7, 10, 11, 9, 8, 13, 14, 12]; 
			return ncon(tensors,connects,con_order=con_order); 
		else 
			error("error: requested environment ($(which_env)) is out of range for current network; please set 'which_env' in range [0-12]") 
		end 


	elseif which_network == 3
		error("error: selected network ($(which_network)) is invalid: network is empty.") 

	elseif which_network == 4
		error("error: selected network ($(which_network)) is invalid: network is empty.") 

	else 
		error("error:requested network is out of range; please set 'which_network' in range [1-4].") 
	end 

end 
