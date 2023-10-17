function MI(rdm)
    
    I_AB = zeros(N,N)
    A,B = size(rdm)
    for a in 1:A
        for b in a+1:B
            rho_AB = copy(rdm[a,b])
            delta_AA = delta(dag(inds(rho_AB)[2]), dag(inds(rho_AB)[1]))
            delta_BB = delta(dag(inds(rho_AB)[4]), dag(inds(rho_AB)[3]))
            @show rho_AB*delta_AA*delta_BB
            rho_B = delta_AA*rho_AB
            rho_A = delta_BB*rho_AB
            
            S_AB = -tr(reshape(Array(rho_AB,inds(rho_AB)[1], inds(rho_AB)[3], inds(rho_AB)[2], inds(rho_AB)[4]),4,4) * log(reshape(Array(rho_AB,inds(rho_AB)[1], inds(rho_AB)[3], inds(rho_AB)[2], inds(rho_AB)[4]),4,4)))
            S_A = -tr(Array(rho_B, inds(rho_B)[1], inds(rho_B)[2]) * log(Array(rho_B, inds(rho_B)[1], inds(rho_B)[2])))
            S_B = -tr(Array(rho_A, inds(rho_A)[1], inds(rho_A)[2]) * log(Array(rho_A, inds(rho_A)[1], inds(rho_A)[2])))
            
            I_AB[a,b] =I_AB[b,a]= S_A + S_B - S_AB
            
        end
    end 
    return I_AB
end

function MI_diag(rdm)
    
    I_AB = zeros(N,N)
    A,B = size(rdm)
    for a in 1:A
        for b in a+1:B
            rho_AB = copy(rdm[a,b]);

            L_inds = (dag(inds(rho_AB)[1]), dag(inds(rho_AB)[3]));
            R_inds = (dag(inds(rho_AB)[2]), dag(inds(rho_AB)[4]));
            delta_AA = delta(dag(inds(rho_AB)[2]), dag(inds(rho_AB)[1]));
            delta_BB = delta(dag(inds(rho_AB)[4]), dag(inds(rho_AB)[3]));

            rho_B = delta_AA*rho_AB;

            rho_A = delta_BB*rho_AB;


            D_AB, _ = eigen(rho_AB, L_inds, R_inds);

            S_AB = 0.0;
            for n=1:dim(D_AB, 1)
                p = D_AB[n,n]
                S_AB -= p * log(p)
            end

            D_A, _ = eigen(rho_A)
            S_A = 0.0
            for n=1:dim(D_A, 1)
                p = D_A[n,n]
                S_A -= p * log(p)
            end

            D_B, _ = eigen(rho_B)
            S_B = 0.0
            for n=1:dim(D_B, 1)
                p = D_B[n,n]
                S_B -= p * log(p)
            end

            println("--------------------------\nRDM\n----------------------")
            @show rho_AB;
            @show rho_A;
            @show rho_B;
            println("---------------------\nCheck norm\n------------------------")
            @show rho_AB*delta_AA*delta_BB;
            println("-----------------------\nEigval\n-----------------------")
            @show D_AB;
            @show D_A;
            @show D_B;

            #S_AB = -tr(reshape(Array(rho_AB,inds(rho_AB)[1], inds(rho_AB)[3], inds(rho_AB)[2], inds(rho_AB)[4]),4,4) * log(reshape(Array(rho_AB,inds(rho_AB)[1], inds(rho_AB)[3], inds(rho_AB)[2], inds(rho_AB)[4]),4,4)))
            #S_A = -tr(Array(rho_B, inds(rho_B)[1], inds(rho_B)[2]) * log(Array(rho_B, inds(rho_B)[1], inds(rho_B)[2])))
            #S_B = -tr(Array(rho_A, inds(rho_A)[1], inds(rho_A)[2]) * log(Array(rho_A, inds(rho_A)[1], inds(rho_A)[2])))
            @show S_A, S_B, S_AB, (a,b);

            I_AB[a,b] =I_AB[b,a]= S_A + S_B - S_AB
            
        end
    end 
    return I_AB
end
