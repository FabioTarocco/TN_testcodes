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
Zero = ITensor([0 1], s)
res = Pup*res
@show res