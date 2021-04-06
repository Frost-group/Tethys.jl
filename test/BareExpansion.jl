@testset "BareExpansion" begin

MAX_ORDER=2
diag=Tethys.Diag(0,1.0,MAX_ORDER,[0.1,0.2,1.0 , 0.3,0.4,1.0], [0,0,0,0,0,0])

H=Tethys.FrohlichHamiltonian(α=1.0, μ=-1.2)

GF=Tethys.BareExpansion(diag, H)
@test GF ≈ -0.7444 atol=0.1 

println("Splish splash: ",diag," ==♓> GF = ", GF)

GF=Tethys.Monte!(diag, H)
println("Monte: ",diag," ==♓>  GF= ", GF)

@time Tethys.Monte!(diag,H)

fout=open("Monte.dat","w")
for i in 1:100
    GF=Tethys.Monte!(diag, H)
    @printf(fout,"%d %f \n",i,GF)
end
close(fout)

end

