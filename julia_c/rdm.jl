
function reduced_rho(psi,i,j)
    #=
        #Move the orthogonality center on the i-th site
        orthogonalize!(psi, i)
    
        #get the dual mps
        psi_dag = dag(psi)
    
        #prime all the links in the dag
        prime(psi_dag, "Link")
    
        #get the link connecting site i-1 and site i
        if i ==1
            rho = psi[i]*prime(psi_dag[i], "Site")
        
        else
            link_i_1 = commonind(psi[i],psi[i-1])
            rho = prime(psi[i],link_i_1) * prime(psi_dag[i], "Site")
            @visualize rho
        end
    
        #all the sites before i-th site are skipped because contract to identity
        #so we can contract dag(psi_i)*psi_i on the link previous site-link
        #during the contraction, it's needed to prime the site i or the index "site" will be contracted between psi and psi_dag
        
    
        #contract all the tensor between the i-th site and the j-th site
        for k in i+1:j-1
            rho *= psi[k]
            rho *= psi_dag[k]
        end
    
        if j == size(psi)[1]
            #if its the last site, the j-th is contracted to rho by link_j-1 leaving Sj
            rho *= psi[j]
            #then the dag site is contracted on link_j-1' and leaving Sj'
            rho *= prime(psi_dag[j], "Site")
        else
            
            #get the link connecting site j and site j+1
            link_j = commonind(psi[j],psi[j+1])
            #rho now has Si, Si', link_j-1, link_j-1' as open indices
            #first we contract the j-th site by link_j-1 leaving open Sj and link_j'
            rho *= prime(psi[j], link_j)
            #the we add the j-th site of the psi_dag by contracting index link_j-1' and link_j' and leaving Sj'
            rho *= prime(psi_dag[j], "Site")
        end
        
        return rho
    =#
     
        #Move the orthogonality center to the i-th site
        orthogonalize!(psi, i);
    
        #create the psi^dag
        psi_dag = dag(psi);
        
        if i==j
            rho = prime(psi_dag[i], "Site") * psi[i];
            @assert 2 == length(inds(rho));
            return rho;
        end
    
        #prime all the virtual index of psi^dag
        prime!(psi_dag, "Link");
    
        #if i is the first site, then we just need to compute the product state of the first site of psi^dag and psi
        if i == 1
            #prime the physical index of psi[i]^dag, otherwise it will be contracted with the physical index of psi[i]
            rho = prime(psi_dag[i],"Site");
            #add the site i to rho
            rho *= psi[i];
        else #if i is not 1, then we need to get link_i-1
            link_i_1 = commonind(psi[i],psi[i-1]);
            #add to rho the i-th site of psi^dag and prime the physical index
            rho = prime(psi_dag[i],"Site");
            #add the i-th site of psi by priming link_i-1
            rho *= prime(psi[i], link_i_1);
        end
        #add to rho all the intermediate sites of psi and psi^dag between site i and j
        for k in i+1:j-1
            rho *= psi_dag[k];
            rho *= psi[k];
        end
    
        #if j is equals to N (last mps site) then we just need to add psi[j]^dag priming the physical index and psi[j]
        if j == size(psi)[1]
            rho *= prime(psi_dag[j], "Site");
            rho *= psi[j];
        else
            #if j is not the last site, we need to get the link index between the j-th site and the j+1-th site
            println("-----------------\nDEGUB\n------------------------");
            @show (j,j+1);
            @show (inds(psi[j]), inds(psi[j+1]));
            link_j = commonind(psi[j], psi[j+1]);
            #add to rho the j-th site of psi^dag priming the physical index and then add the j-th site of psi
            rho *= prime(psi_dag[j], "Site");
            rho *= prime(psi[j], link_j);
        end
    
        @assert 4 == length(inds(rho));
        return rho;
    end



function reduced_rho_matrix(psi)
    el = size(psi)[1]
    density_matrix = Array{Any}(undef, el,el)
    @show size(density_matrix)
    
    for i in 1:el
        for j in i:el
            density_matrix[i,j] = reduced_rho(psi, i, j);
        end
    end
    return density_matrix
end