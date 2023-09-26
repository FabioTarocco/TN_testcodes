using ITensors
using PastaQ


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

i = Index(4, "i")
j = Index(3, "j")
l = Index(4, "l")

A = randomITensor(i,j,l)

trA = A * delta(i,l)
@show trA

T = randomITensor(i,j,l)
U,S,Vh = svd(T,(i,j));

@show norm(U*S*Vh - T)

truncerr = (norm(U*S*Vh - T)/norm(T))^2


i = Index(10, "i")
j = Index(20, "j")
k = Index(20, "k")

T = randomITensor(i,j,k)

U, S, Vh = svd(T, (i,k), cutoff = 1E-2);
@show norm(U*S*Vh - T)
@show (norm(U*S*Vh - T)/norm(T))^2


T = randomITensor(i,j,k)

Q, R = qr(T, (i,k); positive = true);

C = combiner( i,k; tags = "c")

@show C
@show inds(C)

combinedT = C*T
@show inds(combinedT)

combined_index = combinedind(C)

T_ = randomITensor(i,k,i',j)
T_ = C * T_
T_ = noprime(T_)
T_

uncombined_combinedT = dag(C)*combinedT
uncombined_combinedT == T


using ITensors
p0 = QN("P", 0, 2)
p1 = QN("P", 1, 2)
a = QN(("N",0), ("Sz",1,2))
b = QN("Sz", -1)

a = a + b

i = Index(  QN("N", 0) => 1,
            QN("N", 1) => 3, 
            QN("N", 2) =>2; tags = "i")
dag_i = dag(i)


#= 
HARD-CORE BOSON

a|1> = |0>
a^dag|0> = |1>
n|0> = 0|0>
n|1> = 1|1>

(out)|1>, (out)a(in) so it can be contracted (out)a(in)(out)|1> = (out)|0>
=#

s = Index(  QN("N", 0) => 1, 
            QN("N", 1) => 1; tags = "Boson")
a = ITensor(s', dag(s)) 
#quindi 2 stati 0 e 1, con al massimo 1 occupazione, con indice dag(s)(in) per la contrazione e s'(out)
a[s'=>1, s=>2] = 1.0
a_dag = ITensor(s', dag(s))
a_dag[s'=>2, s=>1] = 1.0
n = ITensor(s', dag(s))
n[s'=>2, s=>2] = 1.0



I = Index( QN(0)=>1, QN(1)=>1, tags = "I")
@show I
dim(I)

#indice J con 5 settori con dimensione variabile
J = Index(QN(-2)=>2, QN(-1)=>4, QN(0)=>6, QN(+1)=>4, QN(+2)=>2, tags = "J")

q2 = QN(("N",2),("Sz",2))
QN("N",3) + QN(("Sz",2),("N",2)) == QN(("N",5),("Sz",2))

QN("P",1,2) + QN("P",1,2)

Sz = op("Sz",)


sites = siteind("S=1/2", 2; conserve_qns=false)
