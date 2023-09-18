using ITensors


function heisenberg_mpo(N)
    sites = siteinds("S=1/2", N)
    os = OpSum();
    for i=1:N-1
        os +=       "Sz",i,"Sz",i+1
        os +=   1/2,"S+",i,"S-",i+1
        os +=   1/2,"S-",i,"S+",i+1
    end
    H = MPO(os, sites);
    println(os)
    return H, sites
end
 
N = 10;
H,sites= heisenberg_mpo(N);
println(H)


state = [isodd(n) ? "Up" : "Dn" for n=1:N]
psi0_i = MPS(sites,state);

sweeps = Sweeps(10);
setmaxdim!(sweeps, 10,20,100,200,400,800)
setcutoff!(sweeps,1E-8);
energy, psi0 = dmrg(H, psi0_i, sweeps);


j = 4
A = psi0[j];
@show A


avgSz = expect(psi0, "Sz")
C = correlation_matrix(psi0, "Sz", "Sz")

maxlinkdim(psi0)

psi_truncated = truncate!(psi0; maxdim=2, cutoff = 1E-8)
C_2 = correlation_matrix(psi_truncated, "Sz", "Sz")


a = Index(3, "a")
b = Index(2, "b")
c = Index(4, "c")
d = Index(5, "d")
i = Index(2, "i")
j = Index(6, "j")

A = randomITensor(a,b,c,d)
B = randomITensor(i,d,j)

C = A * B

@show hasinds(C,a,b,c,i,j)


sites = siteinds("S=1/2", N)
os = OpSum();
for j=1:N-1
    os += "Sz",j,"Sz",j+1
    os += 1/2,"S+",j,"S-",j+1
    os += 1/2,"S-",j,"S+",j+1
end
H = MPO(os,sites)
psi_0 = randomMPS(sites; linkdims=10)
sweeps = Sweeps(5)
setmaxdim!(sweeps, 10, 20, 100, 100, 200)
setcutoff!(sweeps, 1E-11)
energy, psi_GS = dmrg(H, psi_0, sweeps)
println("GS eergy = $energy")


i = Index(3, "index_i")
j = Index(2, "index_j")
k = Index(4, "index_k")

T = randomITensor(i,j,k)
@show T

T[i=>1,j=>1,k=>1] = 10.0

T[i=>1,j=>1,k=>1]

M = [1.0 0; 0 1.0]

i = Index(2, "i")
j = Index(2, "j")

A = ITensor(M,i,j)
@show A
a = ITensor([1.0 0], i)

b  = A*a

@show b_ = ITensor([1.0 0], j)
@show b*b_


i = Index(6, "i")
O1 = onehot(i=>6)
@show O1


