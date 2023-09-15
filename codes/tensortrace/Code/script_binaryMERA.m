% script_binaryMERA.m
% ---------------------------------------------------------------------
% Script file demonstrating how the auto-generated code from the
% "TensorTrace" software can be implemented in a tensor network algorithm
% (in this case an optimization of a binary MERA for the ground state of
% the 1D transverse field quantum Ising model on a finite lattice). This
% script calls the "binaryMERA.m" function file as automatically generated
% from the example TensorTrace project "binaryMERA.ttp" included with the
% distribution.
% 
% by Glen Evenbly (www.glenevenbly.com) - last modified 30/8/2019


% set simulation parameters
chi = 6;
chi_p = 4;
n_levels = 3;
n_iterations = 2000;
n_sites = 3*2^(n_levels + 1);

% define hamiltonian
sX = [0,1;1,0];
sY = [0,-1i;1i,0];
sZ = [-1,0;0,1];
ham_s = real(-kron(sX,sX)+0.5*kron(eye(2),sZ)+0.5*kron(eye(2),sZ));
ham_init = 0.5*kron(eye(8),kron(ham_s,eye(2))) + kron(eye(4),kron(ham_s,eye(4))) +...
    0.5*kron(eye(2),kron(ham_s,eye(8)));
chi_b = 4;

% initialize tensors
uDis = cell(n_levels,1); uDis{1} = reshape(eye(chi_b^2,chi_b^2),[chi_b,chi_b,chi_b,chi_b]); 
wIso = cell(n_levels,1); wIso{1} = rand(chi_b,chi_b,chi); 
hamThree = cell(n_levels,1); hamThree{1} = reshape(ham_init - max(eigs(ham_init))*eye(chi_b^3),[chi_b,chi_b,chi_b,chi_b,chi_b,chi_b]);
rhoThree = cell(n_levels,1); rhoThree{1} = rand(chi_b,chi_b,chi_b,chi_b,chi_b,chi_b);
for k = 2:n_levels
    uDis{k} = reshape(eye(chi^2,chi_p^2),[chi,chi,chi_p,chi_p]); 
    wIso{k} = rand(chi_p,chi_p,chi);
    hamThree{k} = rand(chi,chi,chi,chi,chi,chi);
    rhoThree{k} = rand(chi,chi,chi,chi,chi,chi);
end
hamThree{n_levels+1} = rand(chi,chi,chi,chi,chi,chi);
rhoThree{n_levels+1} = rand(chi,chi,chi,chi,chi,chi);

% do optimization iterations
for p = 1:1:n_iterations
    
    % sweep over all levels
    for z = 1:1:n_levels
        if p > 5
            % optimise disentanglers
            tensor_list = {uDis{z}, wIso{z}, hamThree{z}, rhoThree{z+1}};
            uEnv = binaryMERA(tensor_list, 1, 1) + binaryMERA(tensor_list, 1, 2) + binaryMERA(tensor_list, 2, 1) + binaryMERA(tensor_list, 2, 2);
            uSize = size(uEnv);
            [utemp,~,vtemp] = svd(reshape(uEnv,[uSize(1)*uSize(2),uSize(3)*uSize(4)]),0);
            uDis{z} = reshape(-utemp*vtemp',uSize);
        end
        
        % optimise isometries
        tensor_list = {uDis{z}, wIso{z}, hamThree{z}, rhoThree{z+1}};
        Wenv = binaryMERA(tensor_list, 1, 3) + binaryMERA(tensor_list, 1, 4) + binaryMERA(tensor_list, 1, 5) +...
            binaryMERA(tensor_list, 2, 3) + binaryMERA(tensor_list, 2, 4) + binaryMERA(tensor_list, 2, 5);
        wSize = size(Wenv);
        [utemp,~,vtemp] = svd(reshape(Wenv,[wSize(1)*wSize(2),wSize(3)]),0);
        wIso{z} = reshape(-utemp*vtemp',wSize);
        
        % lift Hamiltonian
        tensor_list = {uDis{z}, wIso{z}, hamThree{z}, rhoThree{z+1}};
        hamThree{z+1} = binaryMERA(tensor_list, 1, 12) + binaryMERA(tensor_list, 2, 12);
    end
    
    % diagonalize Hamiltonian
    ham_top = reshape(hamThree{n_levels+1}+permute(hamThree{n_levels+1},[2,3,1,5,6,4])+permute(hamThree{n_levels+1},[3,1,2,6,4,5]),[chi^3,chi^3]);
    [vtemp,dtemp] = eigs(0.5*(ham_top+ham_top'),1,'SA');
    vtemp = vtemp(:)/norm(vtemp);
    rhoThree{n_levels+1} = reshape(vtemp*(vtemp'),[chi,chi,chi,chi,chi,chi]);
    
    % lower the density matrix
    for z = n_levels:-1:1
        tensor_list = {uDis{z}, wIso{z}, hamThree{z}, rhoThree{z+1}};
        rhoThree{z} = 0.5*(binaryMERA(tensor_list, 1, 6) + binaryMERA(tensor_list, 2, 6));
    end
    
    % compute energy and magnetization
    Energy_per_site = trace(reshape(rhoThree{1},[chi_b^3,chi_b^3])*ham_init)/2;
    ExpectX = trace(reshape(rhoThree{1},[chi_b^3,chi_b^3])*kron(eye(32),sX));
    EnExact = (-2/sin(pi/(2*n_sites))) / n_sites;
    EnError = Energy_per_site - EnExact;
    fprintf('Iteration: %3.0d of %3.0d, Energy: %12.12d, Energy Error: %1.4d, XMag: %2.2d\n',p,n_iterations,Energy_per_site,EnError,ExpectX);
end



