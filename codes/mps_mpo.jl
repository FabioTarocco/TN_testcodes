using ITensors

using PastaQ

cutoff = 1E-8
maxddim = 10

i = Index(10, "i")
j = Index(20, "j")
k = Index(20, "k")
l = Index(30, "l")
m = Index(40, "m")
T = randomITensor(i,j,k,l,m)
M = MPS(T, (i,j,k,l,m); cutoff=cutoff, maxdim=maxddim)


d = 2
N = 5
A = randn(d,d,d,d,d)
A = randn(d^N)
sites = siteinds(d, N)
M = MPS(A,sites; cutoff=cutoff, maxdim = maxddim)

N = 10 
s = siteinds(2,N)
@show s
chi = 4
psi = randomMPS(s; linkdims = chi)


N = 4
chi = 4
sites = siteinds("S=1/2", N);
psi = randomMPS(sites; linkdims = chi);
@show psi
magz = expect(psi, "Sz");
for (j, mz) in enumerate(magz)
    println("$j $mz")
end

sum(magz)

W = randomMPO(sites);

ex_W = inner(psi',W,psi);

s = Index(2, "S=1/2")
Z = op("Z",s)

Zp = state("Zp", s)
Zm = state("Zm", s)
@show Z*Zm


using ITensorUnicodePlots
(dag(prime(Zm))*Z*Zm)[]



#Se mi serve lavorare con siti come ad esempio quelli con spin 1/2 allora posso definire
#tipologia specifica di indici

s = Index(2, "S=1/2")
Sz = op("Sz", s)
@show Sz
Z = op("Z", s)
@show Z
q = Index(2, "Qubit")
Z = op("Z", q)

N = 4
sites = siteinds("S=1/2", N)
Sz2 = op("Sz", sites[2])
mps = randomMPS(sites; linkdims=10)
newmps = Sz2 * mps[2]
mps[2] = newmps


s = siteind("S=1/2")
Sz = op("Sz", s)
@show Sz

ITensors.op(::OpName"Pup", ::SiteType"S=1/2") = [1 0;0 0]
Pup = op("Pup",s)
@show Pup
Zero = ITensor([1 0], s)
res = dag(prime(Zero))*Pup*Zero
@show res[]


N = 40
s = siteinds("S=1/2", N)
Pup1 = op("Pup", s[1])
Pup3 = op("Pup", s[3])
@show Pup3


N = 100
sites = siteinds("S=1/2", N)
os = OpSum()
for n=1:N
    os+= "Pup",N
end
P = MPO(os, sites)
mps = randomMPS(sites, linkdims=10)
res = inner(mps', P, mps)

ITensors.state(::StateName"Emp", ::SiteType"Electron") = [1.0, 0, 0, 0]
ITensors.state(::StateName"Up", ::SiteType"Electron") = [0, 1.0, 0, 0]
ITensors.state(::StateName"Dn", ::SiteType"Electron") = [0, 0, 1.0, 0]
ITensors.state(::StateName"UpDn", ::SiteType"Electron") = [0, 0, 0, 1.0]
ITensors.state(::StateName"+", ::SiteType"Electron") = [0, 1/sqrt(2), 1/sqrt(2), 0]

s = siteind("Electron")
plus = state("+", s)
@show tr(plus'*plus)


#SISTEMI SPIN 3/2
ITensors.space(::SiteType"S=3/2") = 4
ITensors.op(::OpName"Sz", ::SiteType"S=3/2") = 
    [+3/2 0 0 0
    0 +1/2 0 0
    0 0 -1/2 0
    0 0 0 -3/2]
ITensors.op(::OpName"S+", ::SiteType"S=3/2") = 
    [0 sqrt(3) 0 0
    0 0 2 0
    0 0 sqrt(2) 0
    0 0 0 0]
ITensors.op(::OpName"S-", ::SiteType"S=3/2") = 
    [0 0 0 0
    sqrt(3) 0 0 0
    0 2 0 0
    0 0 sqrt(3) 0]

s = siteind("S=3/2")
N = 100
sites = siteinds("S=3/2", N)

s = Index(4, "S=3/2")
Sz = op("Sz", s)
@show Sz
println(Sz)
Sz1 = op("Sz", sites[1])
Sz3 = op("S+", sites[3])


N = 100
sites = siteinds("S=1", N)
#Spin 1/2 -> dim=2, Spin 1 -> dim=3, spin 3/2 -> dim = 3/2
os = OpSum()
for j=1:N-1
    os += "Sz",j,"Sz",j+1
    os += 1/2, "S+",j,"S-",j+1
    os += 1/2, "S-",j,"S+",j+1
end
H = MPO(os, sites);
psi_0 = randomMPS(sites,linkdims=10);
nsweeps = 5
maxdim =[10,20,100,100,200]
cutoff = [1E-10]
using Profile
en, psi_gs = dmrg(H,psi_0; nsweeps, maxdim, cutoff)
@show norm(psi_gs)
@show sites[1]

N = 20
sites = siteinds("S=1", N; conserve_qns = true)
initial_conf = [isodd(n) ? "Up" : "Dn" for n=1:N]
psi_0 = randomMPS(sites, initial_conf)
@show flux(psi_0)
initial_all_up = ["Up" for n=1:N]
psi_0_all = randomMPS(sites, initial_all_up)
@show flux(psi_0_all)

N = 100
sites = siteinds("S=1",N; conserve_qns=true)
os = OpSum()
for j=1:N-1
    os += "Sz",j,"Sz",j+1
    os += 1/2, "S+",j,"S-",j+1
    os += 1/2, "S-",j,"S+",j+1
end
H = MPO(os, sites)
initial_updn = [isodd(n) ? "Up" : "Dn" for n=1:N]
psi_0 = randomMPS(sites, initial_updn)
@show flux(psi_0)

nsweeps = 6
maxdim = [10,20,100,100,200]
cutoff = [1E-10]

en, psi_dmrg = dmrg(H,psi_0; nsweeps, maxdim, cutoff)

#SUZUKI-TROTTER DECOMPOSITION for TIME EVOLUTION BLOCK DECIMATION (TEBD)
N = 100
cutoff =1E-8
tau = 0.1 #time rate of every single step
ttotal = 5.0 #total time for the evolution

s = siteinds("S=1/2", N; conserve_qns = true)

gates = ITensor[]
for j in 1:(N-1)
    s1 = s[j]
    s2 = s[j+1]
    hj = 
        op("Sz",s1)*op("Sz",s2) +
        1/2 * op("S+", s1)*op("S-",s2)+        
        1/2 * op("S-", s1)*op("S+",s2)
        #definizione del sigolo blocco locale dell'evoluzione

        Gj = exp(-im * tau / 2 * hj)
        push!(gates, Gj)
end

append!(gates, reverse(gates))
initial_conf = [isodd(n) ? "Up" : "Dn" for n=1:N]
psi = MPS(s,initial_conf)
c = div(N,2)

for t in 0.0:tau:ttotal
    Sz = expect(psi, "Sz"; sites = c)
    println("$t $Sz")
    t≈ttotal && break
    psi = apply(gates, psi; cutoff)
    normalize!(psi)
end



###
#----------------------DIFFERENT SYSTEM MPO DMRG------------------------
###

#---------------------HEISENBERG 1D  (S=1)      ------------------------
using Printf
using Random

Random.seed!(1234)
N = 100
sites = siteinds("S=1", N)
os = OpSum()
for j in 1:(N-1)
    os += "Sz",j,"Sz",j+1
    os += 1/2,"S+",j,"S-",j+1
    os += 1/2,"S-",j,"S+",j+1
end

H = MPO(os, sites)
psi_0 = randomMPS(sites, 10)
nsweeps = 5
maxdim = [10,20,100,100,200]
cutoff = [1E-11]
en_dmrg, psi_dmrg = dmrg(H, psi_0; nsweeps, maxdim, cutoff)
@printf("Final energy = %1.12f\n", en_dmrg)


#-------------------HEISENBERG 1D (S=1/2) ---------------------
N = 100
sites = siteinds("S=1/2", N)
os = OpSum()
for j in 1:(N-1)
    os += "Sz",j,"Sz",j+1
    os += 1/2,"S+",j,"S-",j+1
    os += 1/2,"S-",j,"S+",j+1
end

H = MPO(os, sites)
psi_0 = randomMPS(sites, 10)
nsweeps = 5
maxdim = [10,20,100,100,200]
cutoff = [1E-11]
en_dmrg, psi_dmrg = dmrg(H, psi_0; nsweeps, maxdim, cutoff)
@printf("Final energy = %1.12f\n", en_dmrg)



#---------------------HEISENBERG 1D  (S=1)  spin preserving ------------------------

N = 100
sites = siteinds("S=1", N; conserve_qns=true)

os = OpSum()
for j in 1:(N-1)
    os += "Sz",j,"Sz",j+1
    os += 1/2,"S+",j,"S-",j+1
    os += 1/2,"S-",j,"S+",j+1
end
H = MPO(os, sites)
initial_conf = [isodd(n) ? "Up" : "Dn" for n in 1:N]
psi_0 = randomMPS(sites, initial_conf, 10)

nsweeps = 5
maxdim = [10, 20, 100, 100, 200]
cutoff = 1E-11

en_dmrg, psi_dmrg = dmrg(H, psi_0; nsweeps, maxdim, cutoff)
@printf("Final energy = %1.12f\n", en_dmrg)


#---------------------EXTENDED HUBBARD 1D  (S=1/2)  spin preserving ------------------------


N = 20
Npart = 10
t1 = 1.0
t2 = 0.2
U = 1.0
V1 = 0.5

sites = siteinds("Electron", N; conserve_qns = true)
os = OpSum()
for b in 1:(N-1)
    os += -t1,  "Cdagup",b,"Cup",b+1
    os += -t1,  "Cup",b,"Cdagup",b+1
    os += -t1,  "Cdagdn",b,"Cdn",b+1
    os += -t1,  "Cdn",b,"Cdagdn",b+1
    os += V1,   "Ntot",b,"Ntot",b+1
end
for b in 1:(N-2)
    os += -t2,  "Cdagup",b,"Cup",b+2
    os += -t2,  "Cup",b,"Cdagup",b+2
    os += -t2,  "Cdagdn",b,"Cdn",b+2
    os += -t2,  "Cdn",b,"Cdagdn",b+2
end
for i in 1:N
    os += U, "Nupdn",i
end
H = MPO(os, sites)

nsweeps = 6
maxdim = [50,100,200,400,800,800]
cutoff = [1E-12]

initial_conf = ["Emp" for n in 1:N]
p = Npart
for i in N:-1:1
    if p>i
        println("Doubly occupied site $i")
        initial_conf[i] = "UpDn"
        p -= 2
    elseif p > 0
        println("Singly occupied site $i")
        initial_conf[i] = (isodd(i) ? "Up" : "Dn")
        p -= 1
    end
end

@show initial_conf

psi_0 = randomMPS(sites, initial_conf, 10)
@show flux(psi_0)

en_dmrg, psi_dmrg = dmrg(H, psi_0; nsweeps, maxdim, cutoff)


#---------------------HEISENBERG 2D  (S=1/2)  spin preserving ------------------------

# definiamo prima il reticolo di interesse (lattice)
# il quale ci permette di ottenere i corretti termini di coupling

Ny = 6 
Nx = 12
N = Nx * Ny #numero totale di siti

sites = siteinds("S=1/2", N; conserve_qns = true)
lattice = square_lattice(Nx, Ny; yperiodic = false) #yperiodic per evitare al formazione di un cilindro
# in ogni caso il lattice numera in questo modo     4 - 5 - 6  e via così
#                                                   |
#                                                   +-------+
#                                                           |
#                                                   1 - 2 - 3

@show lattice

os = OpSum()
for b in lattice 
    os += "Sx", b.s1, "Sx", b.s2
    os += "Sy", b.s1, "Sy", b.s2 
    os += "Sz", b.s1, "Sz", b.s2
end
H_spin = MPO(os, sites)

os_1 = OpSum()
for b in lattice 
    os_1 += 1/2, "S+", b.s1, "S-", b.s2
    os_1 += 1/2, "S-", b.s1, "S+", b.s2 
    os_1 += "Sz", b.s1, "Sz", b.s2
end
H_pm = MPO(os_1, sites)

initial_condition = [isodd(n) ? "Up" : "Dn" for n in 1:N]
psi_0 = randomMPS(sites, initial_condition, 20)

nsweeps = 10
maxdim = [20, 60, 100, 100, 200, 400, 800]
cutoff = [1E-8]


en_dmrg, psi_dmrg = dmrg(H, psi_dmrg; nsweeps, maxdim, cutoff)

#------------------------------- 2D HUBBARD conserve num particles -------------------

Nx = 6
Ny = 3
U = 4.0
t = 1.0
N = Nx * Ny
nsweeps = 10
maxim = [100,200,400,800,1600]
cutoff = [1E-6]
sites = siteinds("Electron", N; conserve_qns = true)
lattice = square_lattice(Nx, Ny; yperdioc = true)
@show lattice[2]

os = OpSum()
for b in lattice
    os += -t, "Cdagup", b.s1, "Cup", b.s2
    os += -t, "Cdagup", b.s2, "Cup", b.s1
    os += -t, "Cdagdn", b.s1, "Cdn", b.s2
    os += -t, "Cdagdn", b.s2, "Cdn", b.s1
end
for i in 1:N
    os += U, "Nupdn", i
end
H = MPO(os, sites)
initial_conf = [isodd(n) ? "Up" : "Dn" for n in 1:N]

psi_random = randomMPS(sites, initial_conf)
energy, psi_dmrg = dmrg(H, psi_random; nsweeps, maxim, cutoff)
@show t,U
@show flux(psi_dmrg)
@show maxlinkdim(psi_dmrg)
@show energy 

