function tensor_out = binaryMERA(tensors_in, which_network, which_env) 
% function tensor_out = binaryMERA({u, w, ham, rho}, which_network, which_env) 
% 
% Auto-generated network contraction function from 'TensorTrace' software, 
% see (www.tensortrace.com) for details, (c) Glen Evenbly, 2019. 
% Requires the network contraction routine 'ncon' to be in the working directory, 
% (included in the TensorTrace install directory, or also found at 'www.tensortrace.com'). 
% 
% Input variables: 
% ------------------------------- 
% 1st input 'tensors_in' is a cell array of the unique tensors: tensors_in = {u, w, ham, rho}; 
% 2nd input 'which_network' dictates which network to evaluate (1-4) 
% 3rd input 'which_env' allows one to specify a tensor environment to evaluate from a closed tensor network. 
% 	 -set 'which_env = 0' to evaluate the scalar from a closed network (i.e. no environment). 
% 	 -set 'which_env = n' to evaluate environment of the nth tensor from the closed network. 
% 
% General project info: 
% ------------------------------- 
% Generated on: 30/6/2020
% Generated from: C:\Users\gevenbly3\Desktop\TensorTrace\SaveData\
% Index dims: chi = 8,chi_p = 6

% --- Info for Network-1 --- 
% -------------------------- 
% network is CLOSED 
% total contraction cost (in scalar multiplications): 8.86*10^7
% contraction order: (((((T1*((T2*(T6*T8))*(T7*T10)))*T4)*(T5*T11))*(T3*T9))*T12)
% 1st leading order cost: (chi^4)*(chi_p^5)
% 2nd leading order cost: (chi^3)*(chi_p^6)
% tensors_N1 = {u,u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w),rho}; 
% connects_N1 = {[1,3,10,11],[4,7,12,13],[8,10,17],[11,12,18],[13,14,19],[2,5,6,3,4,7],[1,2,9,23],[5,6,16,15],[8,9,22],[23,16,21],[15,14,20],[22,21,20,17,18,19]}; 
% dims_N1 = {[chi,chi,chi_p,chi_p],[chi,chi,chi_p,chi_p],[chi_p,chi_p,chi],[chi_p,chi_p,chi],[chi_p,chi_p,chi],[chi,chi,chi,chi,chi,chi],[chi,chi,chi_p,chi_p],[chi,chi,chi_p,chi_p],[chi_p,chi_p,chi],[chi_p,chi_p,chi],[chi_p,chi_p,chi],[chi,chi,chi,chi,chi,chi]}; 
% con_order_N1 = [5, 6, 14, 8, 23, 4, 7, 2, 16, 1, 3, 11, 12, 13, 15, 10, 9, 21, 18, 19, 20, 17, 22]; 

% --- Info for Network-2 --- 
% -------------------------- 
% network is CLOSED 
% total contraction cost (in scalar multiplications): 8.86*10^7
% contraction order: (((((((T1*(T6*T7))*(T2*T4))*T8)*T10)*(T3*T9))*(T5*T11))*T12)
% 1st leading order cost: (chi^4)*(chi_p^5)
% 2nd leading order cost: (chi^3)*(chi_p^6)
% tensors_N2 = {u,u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w),rho}; 
% connects_N2 = {[2,3,9,10],[6,23,11,12],[7,9,16],[10,11,17],[12,13,18],[1,4,5,2,3,6],[1,4,8,22],[5,23,15,14],[7,8,21],[22,15,20],[14,13,19],[21,20,19,16,17,18]}; 
% dims_N2 = {[chi,chi,chi_p,chi_p],[chi,chi,chi_p,chi_p],[chi_p,chi_p,chi],[chi_p,chi_p,chi],[chi_p,chi_p,chi],[chi,chi,chi,chi,chi,chi],[chi,chi,chi_p,chi_p],[chi,chi,chi_p,chi_p],[chi_p,chi_p,chi],[chi_p,chi_p,chi],[chi_p,chi_p,chi],[chi,chi,chi,chi,chi,chi]}; 
% con_order_N2 = [1, 4, 11, 2, 3, 10, 6, 5, 23, 13, 7, 22, 15, 9, 8, 12, 14, 17, 20, 16, 21, 18, 19]; 

% --- Info for Network-3 --- 
% -------------------------- 
% -- Network is empty. -- 

% --- Info for Network-4 --- 
% -------------------------- 
% -- Network is empty. -- 


% Contraction Code: 
% ------------------------------- 

if isempty(tensors_in)% auto-generate random tensors of default dims 
	chi = 8; 
	chi_p = 6; 
	u = rand(chi,chi,chi_p,chi_p); 
	w = rand(chi_p,chi_p,chi); 
	ham = rand(chi,chi,chi,chi,chi,chi); 
	rho = rand(chi,chi,chi,chi,chi,chi); 
else
	u = tensors_in{1}; 
	w = tensors_in{2}; 
	ham = tensors_in{3}; 
	rho = tensors_in{4}; 
end 

if ~exist('which_env','var') 
	which_env = 0; 
end 

if which_network == 1
	if which_env == 0 
		% TTv1.0.5.0M$;@<EHE<?H?@TE`ET?`?B09<963BH9T9N3B`9l9f3>HQTQ`QHKTK`K@<WHW<]H]@TW`WT]`]B0c<c6iBHcTcNiB`clcfiF6'N'f'6-N-f-')++))+++''+''''++++(((LV,HQ^LO=]'IQ^:HUYqHQ^VJ$; 
		tensors = {u,u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w),rho}; 
		connects = {[1,3,10,11],[4,7,12,13],[8,10,17],[11,12,18],[13,14,19],[2,5,6,3,4,7],[1,2,9,23],[5,6,16,15],[8,9,22],[23,16,21],[15,14,20],[22,21,20,17,18,19]}; 
		con_order = [5, 6, 14, 8, 23, 4, 7, 2, 16, 1, 3, 11, 12, 13, 15, 10, 9, 21, 18, 19, 20, 17, 22]; 
		tensor_out = ncon(tensors,connects,con_order); 
	elseif which_env == 1 
		tensors = {u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w),rho}; 
		connects = {[4,7,12,13],[8,-3,17],[-4,12,18],[13,14,19],[2,5,6,-2,4,7],[-1,2,9,23],[5,6,16,15],[8,9,22],[23,16,21],[15,14,20],[22,21,20,17,18,19]}; 
		con_order = [5, 6, 14, 8, 23, 4, 7, 2, 16, 17, 22, 20, 19, 18, 12, 13, 15, 9, 21]; 
		tensor_out = ncon(tensors,connects,con_order); 
	elseif which_env == 2 
		tensors = {u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w),rho}; 
		connects = {[1,3,10,11],[8,10,17],[11,-3,18],[-4,14,19],[2,5,6,3,-1,-2],[1,2,9,23],[5,6,16,15],[8,9,22],[23,16,21],[15,14,20],[22,21,20,17,18,19]}; 
		con_order = [5, 6, 14, 8, 23, 17, 22, 20, 19, 18, 10, 11, 1, 9, 21, 3, 15, 2, 16]; 
		tensor_out = ncon(tensors,connects,con_order); 
	elseif which_env == 3 
		tensors = {u,u,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w),rho}; 
		connects = {[1,3,-2,11],[4,7,12,13],[11,12,18],[13,14,19],[2,5,6,3,4,7],[1,2,9,23],[5,6,16,15],[-1,9,22],[23,16,21],[15,14,20],[22,21,20,-3,18,19]}; 
		con_order = [5, 6, 14, 23, 4, 7, 2, 16, 1, 3, 11, 12, 13, 15, 21, 18, 19, 20, 9, 22]; 
		tensor_out = ncon(tensors,connects,con_order); 
	elseif which_env == 4 
		tensors = {u,u,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w),rho}; 
		connects = {[1,3,10,-1],[4,7,-2,13],[8,10,17],[13,14,19],[2,5,6,3,4,7],[1,2,9,23],[5,6,16,15],[8,9,22],[23,16,21],[15,14,20],[22,21,20,17,-3,19]}; 
		con_order = [5, 6, 14, 8, 23, 4, 7, 2, 16, 1, 3, 17, 22, 20, 19, 10, 13, 15, 9, 21]; 
		tensor_out = ncon(tensors,connects,con_order); 
	elseif which_env == 5 
		tensors = {u,u,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w),rho}; 
		connects = {[1,3,10,11],[4,7,12,-1],[8,10,17],[11,12,18],[2,5,6,3,4,7],[1,2,9,23],[5,6,16,15],[8,9,22],[23,16,21],[15,-2,20],[22,21,20,17,18,-3]}; 
		con_order = [5, 6, 8, 23, 4, 7, 2, 16, 1, 3, 11, 12, 17, 22, 10, 9, 21, 18, 15, 20]; 
		tensor_out = ncon(tensors,connects,con_order); 
	elseif which_env == 6 
		tensors = {u,u,w,w,w,conj(u),conj(u),conj(w),conj(w),conj(w),rho}; 
		connects = {[1,-4,10,11],[-5,-6,12,13],[8,10,17],[11,12,18],[13,14,19],[1,-1,9,23],[-2,-3,16,15],[8,9,22],[23,16,21],[15,14,20],[22,21,20,17,18,19]}; 
		con_order = [14, 8, 23, 17, 22, 20, 19, 18, 10, 11, 1, 9, 21, 13, 12, 15, 16]; 
		tensor_out = ncon(tensors,connects,con_order); 
	elseif which_env == 7 
		tensors = {u,u,w,w,w,ham,conj(u),conj(w),conj(w),conj(w),rho}; 
		connects = {[-1,3,10,11],[4,7,12,13],[8,10,17],[11,12,18],[13,14,19],[-2,5,6,3,4,7],[5,6,16,15],[8,-3,22],[-4,16,21],[15,14,20],[22,21,20,17,18,19]}; 
		con_order = [5, 6, 14, 8, 4, 7, 17, 22, 20, 19, 18, 10, 11, 3, 13, 15, 12, 21, 16]; 
		tensor_out = ncon(tensors,connects,con_order); 
	elseif which_env == 8 
		tensors = {u,u,w,w,w,ham,conj(u),conj(w),conj(w),conj(w),rho}; 
		connects = {[1,3,10,11],[4,7,12,13],[8,10,17],[11,12,18],[13,14,19],[2,-1,-2,3,4,7],[1,2,9,23],[8,9,22],[23,-3,21],[-4,14,20],[22,21,20,17,18,19]}; 
		con_order = [14, 8, 23, 17, 22, 20, 19, 18, 10, 11, 1, 9, 21, 13, 12, 3, 2, 4, 7]; 
		tensor_out = ncon(tensors,connects,con_order); 
	elseif which_env == 9 
		tensors = {u,u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),rho}; 
		connects = {[1,3,10,11],[4,7,12,13],[-1,10,17],[11,12,18],[13,14,19],[2,5,6,3,4,7],[1,2,-2,23],[5,6,16,15],[23,16,21],[15,14,20],[-3,21,20,17,18,19]}; 
		con_order = [5, 6, 14, 23, 4, 7, 2, 16, 1, 3, 11, 12, 13, 15, 21, 18, 19, 20, 10, 17]; 
		tensor_out = ncon(tensors,connects,con_order); 
	elseif which_env == 10 
		tensors = {u,u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),rho}; 
		connects = {[1,3,10,11],[4,7,12,13],[8,10,17],[11,12,18],[13,14,19],[2,5,6,3,4,7],[1,2,9,-1],[5,6,-2,15],[8,9,22],[15,14,20],[22,-3,20,17,18,19]}; 
		con_order = [5, 6, 14, 8, 4, 7, 17, 22, 20, 19, 18, 10, 11, 3, 13, 15, 12, 1, 9, 2]; 
		tensor_out = ncon(tensors,connects,con_order); 
	elseif which_env == 11 
		tensors = {u,u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),rho}; 
		connects = {[1,3,10,11],[4,7,12,13],[8,10,17],[11,12,18],[13,-2,19],[2,5,6,3,4,7],[1,2,9,23],[5,6,16,-1],[8,9,22],[23,16,21],[22,21,-3,17,18,19]}; 
		con_order = [5, 6, 8, 23, 4, 7, 2, 16, 1, 3, 11, 12, 17, 22, 10, 9, 21, 18, 13, 19]; 
		tensor_out = ncon(tensors,connects,con_order); 
	elseif which_env == 12 
		tensors = {u,u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w)}; 
		connects = {[1,3,10,11],[4,7,12,13],[8,10,-4],[11,12,-5],[13,14,-6],[2,5,6,3,4,7],[1,2,9,23],[5,6,16,15],[8,9,-1],[23,16,-2],[15,14,-3]}; 
		con_order = [5, 6, 14, 8, 23, 4, 7, 2, 16, 1, 3, 11, 12, 13, 15, 10, 9]; 
		tensor_out = ncon(tensors,connects,con_order); 
	else 
		error("error: requested environment (%d) is out of range for current network; please set 'which_env' in range [0-12]",which_env); 
	end 


elseif which_network == 2
	if which_env == 0 
		% TTv1.0.5.0M$;@<EHE<?H?@TE`ET?`?B09<963BH9T9N3B`9l9f3><QHQTQ<KHKTK@<WHW<]H]@TW`WT]`]B0c<c6iBHcTcNiB`clcfiF6'N'f'6-N-f-))++)'+++''+''''++++(((LV,HQ^LO=]'IQ^:HUYqHQ^VJ$; 
		tensors = {u,u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w),rho}; 
		connects = {[2,3,9,10],[6,23,11,12],[7,9,16],[10,11,17],[12,13,18],[1,4,5,2,3,6],[1,4,8,22],[5,23,15,14],[7,8,21],[22,15,20],[14,13,19],[21,20,19,16,17,18]}; 
		con_order = [1, 4, 11, 2, 3, 10, 6, 5, 23, 13, 7, 22, 15, 9, 8, 12, 14, 17, 20, 16, 21, 18, 19]; 
		tensor_out = ncon(tensors,connects,con_order); 
	elseif which_env == 1 
		tensors = {u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w),rho}; 
		connects = {[6,23,11,12],[7,-3,16],[-4,11,17],[12,13,18],[1,4,5,-1,-2,6],[1,4,8,22],[5,23,15,14],[7,8,21],[22,15,20],[14,13,19],[21,20,19,16,17,18]}; 
		con_order = [1, 4, 11, 13, 7, 18, 19, 16, 21, 20, 14, 15, 23, 12, 17, 6, 8, 22, 5]; 
		tensor_out = ncon(tensors,connects,con_order); 
	elseif which_env == 2 
		tensors = {u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w),rho}; 
		connects = {[2,3,9,10],[7,9,16],[10,-3,17],[-4,13,18],[1,4,5,2,3,-1],[1,4,8,22],[5,-2,15,14],[7,8,21],[22,15,20],[14,13,19],[21,20,19,16,17,18]}; 
		con_order = [1, 4, 2, 3, 13, 7, 18, 19, 16, 21, 20, 14, 15, 9, 5, 8, 22, 10, 17]; 
		tensor_out = ncon(tensors,connects,con_order); 
	elseif which_env == 3 
		tensors = {u,u,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w),rho}; 
		connects = {[2,3,-2,10],[6,23,11,12],[10,11,17],[12,13,18],[1,4,5,2,3,6],[1,4,8,22],[5,23,15,14],[-1,8,21],[22,15,20],[14,13,19],[21,20,19,-3,17,18]}; 
		con_order = [1, 4, 11, 2, 3, 10, 6, 5, 23, 13, 22, 15, 18, 19, 12, 17, 14, 20, 8, 21]; 
		tensor_out = ncon(tensors,connects,con_order); 
	elseif which_env == 4 
		tensors = {u,u,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w),rho}; 
		connects = {[2,3,9,-1],[6,23,-2,12],[7,9,16],[12,13,18],[1,4,5,2,3,6],[1,4,8,22],[5,23,15,14],[7,8,21],[22,15,20],[14,13,19],[21,20,19,16,-3,18]}; 
		con_order = [1, 4, 2, 3, 13, 7, 18, 19, 16, 21, 20, 14, 15, 9, 5, 8, 22, 6, 12, 23]; 
		tensor_out = ncon(tensors,connects,con_order); 
	elseif which_env == 5 
		tensors = {u,u,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w),rho}; 
		connects = {[2,3,9,10],[6,23,11,-1],[7,9,16],[10,11,17],[1,4,5,2,3,6],[1,4,8,22],[5,23,15,14],[7,8,21],[22,15,20],[14,-2,19],[21,20,19,16,17,-3]}; 
		con_order = [1, 4, 11, 2, 3, 10, 6, 5, 23, 7, 22, 15, 9, 8, 17, 20, 16, 21, 14, 19]; 
		tensor_out = ncon(tensors,connects,con_order); 
	elseif which_env == 6 
		tensors = {u,u,w,w,w,conj(u),conj(u),conj(w),conj(w),conj(w),rho}; 
		connects = {[-4,-5,9,10],[-6,23,11,12],[7,9,16],[10,11,17],[12,13,18],[-1,-2,8,22],[-3,23,15,14],[7,8,21],[22,15,20],[14,13,19],[21,20,19,16,17,18]}; 
		con_order = [11, 13, 7, 18, 19, 16, 21, 20, 14, 15, 23, 12, 17, 9, 10, 8, 22]; 
		tensor_out = ncon(tensors,connects,con_order); 
	elseif which_env == 7 
		tensors = {u,u,w,w,w,ham,conj(u),conj(w),conj(w),conj(w),rho}; 
		connects = {[2,3,9,10],[6,23,11,12],[7,9,16],[10,11,17],[12,13,18],[-1,-2,5,2,3,6],[5,23,15,14],[7,-3,21],[-4,15,20],[14,13,19],[21,20,19,16,17,18]}; 
		con_order = [11, 13, 7, 18, 19, 16, 21, 20, 14, 15, 23, 12, 17, 9, 10, 2, 3, 6, 5]; 
		tensor_out = ncon(tensors,connects,con_order); 
	elseif which_env == 8 
		tensors = {u,u,w,w,w,ham,conj(u),conj(w),conj(w),conj(w),rho}; 
		connects = {[2,3,9,10],[6,-2,11,12],[7,9,16],[10,11,17],[12,13,18],[1,4,-1,2,3,6],[1,4,8,22],[7,8,21],[22,-3,20],[-4,13,19],[21,20,19,16,17,18]}; 
		con_order = [1, 4, 11, 2, 3, 10, 6, 13, 7, 18, 19, 16, 21, 20, 9, 8, 22, 12, 17]; 
		tensor_out = ncon(tensors,connects,con_order); 
	elseif which_env == 9 
		tensors = {u,u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),rho}; 
		connects = {[2,3,9,10],[6,23,11,12],[-1,9,16],[10,11,17],[12,13,18],[1,4,5,2,3,6],[1,4,-2,22],[5,23,15,14],[22,15,20],[14,13,19],[-3,20,19,16,17,18]}; 
		con_order = [1, 4, 11, 2, 3, 10, 6, 5, 23, 13, 22, 15, 18, 19, 12, 17, 14, 20, 9, 16]; 
		tensor_out = ncon(tensors,connects,con_order); 
	elseif which_env == 10 
		tensors = {u,u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),rho}; 
		connects = {[2,3,9,10],[6,23,11,12],[7,9,16],[10,11,17],[12,13,18],[1,4,5,2,3,6],[1,4,8,-1],[5,23,-2,14],[7,8,21],[14,13,19],[21,-3,19,16,17,18]}; 
		con_order = [1, 4, 11, 2, 3, 10, 6, 5, 23, 13, 7, 18, 19, 16, 21, 9, 8, 12, 17, 14]; 
		tensor_out = ncon(tensors,connects,con_order); 
	elseif which_env == 11 
		tensors = {u,u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),rho}; 
		connects = {[2,3,9,10],[6,23,11,12],[7,9,16],[10,11,17],[12,-2,18],[1,4,5,2,3,6],[1,4,8,22],[5,23,15,-1],[7,8,21],[22,15,20],[21,20,-3,16,17,18]}; 
		con_order = [1, 4, 11, 2, 3, 10, 6, 5, 23, 7, 22, 15, 9, 8, 17, 20, 16, 21, 12, 18]; 
		tensor_out = ncon(tensors,connects,con_order); 
	elseif which_env == 12 
		tensors = {u,u,w,w,w,ham,conj(u),conj(u),conj(w),conj(w),conj(w)}; 
		connects = {[2,3,9,10],[6,23,11,12],[7,9,-4],[10,11,-5],[12,13,-6],[1,4,5,2,3,6],[1,4,8,22],[5,23,15,14],[7,8,-1],[22,15,-2],[14,13,-3]}; 
		con_order = [1, 4, 11, 2, 3, 10, 6, 5, 23, 13, 7, 22, 15, 9, 8, 12, 14]; 
		tensor_out = ncon(tensors,connects,con_order); 
	else 
		error("error: requested environment (%d) is out of range for current network; please set 'which_env' in range [0-12]",which_env); 
	end 


elseif which_network == 3
	error("error: selected network (%d) is invalid: network is empty.",which_network); 

elseif which_network == 4
	error("error: selected network (%d) is invalid: network is empty.",which_network); 

else 
	error("error: requested network is out of range; please set 'which_network' in range [1-4]."); 
end 

