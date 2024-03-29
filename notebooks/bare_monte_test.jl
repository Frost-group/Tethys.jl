### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ 311bed50-ee23-4546-9140-a8b67e76d755
# setup package environment over Tethys packages
begin
    import Pkg
    Pkg.activate(Base.current_project())

    using Tethys
	
	using Plots, StatsPlots
end


# ╔═╡ 89a4e2a6-d7eb-11ea-3b22-fb24e94f2aa9
using BenchmarkTools

# ╔═╡ da1bb8e9-92b0-4c1d-a48d-fa1046f441c7
gr(format=:png, size=(800,600)) # Force raster / PNG, as otherwise millions of data points get sent to web browser ^_^

# ╔═╡ 141e94f8-7362-11eb-3606-6360ade0fefc
Ĥ=Tethys.FrohlichHamiltonian(α=5.0, μ=-6)

# ╔═╡ 9a8c77a0-a593-11ea-09e0-5177021a936e
begin
	MAX_ORDER=2
	diag1=Tethys.Diag(0,1.0,MAX_ORDER,[0.1,0.2,1.0 , 0.3,0.4,1.0], [0.0,0.0,0.0,0.0,0.0,0.0])
	GF1=Tethys.BareExpansion(diag1, Ĥ)
end

# ╔═╡ 6c3748a8-a5a7-11ea-2490-d5e6dce272a8
diag1

# ╔═╡ ddcd5f98-a593-11ea-1b0f-6f6067269fc4
GF2=Tethys.Monte!(diag1, Ĥ)

# ╔═╡ 69b81526-a5a7-11ea-272f-21121dbbfbb9
diag1

# ╔═╡ ed0ad5da-a593-11ea-0d19-1fa974170c9d
begin
	diag=Tethys.Diag(0,1.0,MAX_ORDER,[0.1,0.2,1.0 , 0.3,0.4,1.0], [0,0,0,0,0,0])
	MC=[]
	for i in 1:1_000_000
		GF=Tethys.Monte!(diag, Ĥ)
		append!(MC,GF)
	end
end

# ╔═╡ 6bef50b2-a594-11ea-3ab0-dd19b06da250
MC

# ╔═╡ 6f29880e-a594-11ea-3a1e-b18693818f43
plot(MC)

# ╔═╡ e8d22018-a593-11ea-2d09-c921536e9287
density(MC)

# ╔═╡ d08f9624-a599-11ea-3af6-5f6dc915480f
diag

# ╔═╡ d16e91d0-a599-11ea-23f1-1963e0d078bb
Tethys.BareExpansion(diag, Ĥ)

# ╔═╡ abd6ba19-16c4-4844-98bc-b1c7ee5716bd
# Benchmarking below here...

# ╔═╡ a93ba57a-d7eb-11ea-26c8-1180ccb530d5
@benchmark Tethys.Monte!(diag1, H)

# ╔═╡ 4bd42778-d7ec-11ea-095b-8f80fca95014
@benchmark Tethys.BareExpansion(diag, H)

# ╔═╡ c2e29b42-7362-11eb-10e1-e1653cdd521e
@benchmark(exp(2.0))

# ╔═╡ Cell order:
# ╠═311bed50-ee23-4546-9140-a8b67e76d755
# ╠═da1bb8e9-92b0-4c1d-a48d-fa1046f441c7
# ╠═141e94f8-7362-11eb-3606-6360ade0fefc
# ╠═9a8c77a0-a593-11ea-09e0-5177021a936e
# ╠═6c3748a8-a5a7-11ea-2490-d5e6dce272a8
# ╠═ddcd5f98-a593-11ea-1b0f-6f6067269fc4
# ╠═69b81526-a5a7-11ea-272f-21121dbbfbb9
# ╠═ed0ad5da-a593-11ea-0d19-1fa974170c9d
# ╠═6bef50b2-a594-11ea-3ab0-dd19b06da250
# ╠═6f29880e-a594-11ea-3a1e-b18693818f43
# ╠═e8d22018-a593-11ea-2d09-c921536e9287
# ╠═d08f9624-a599-11ea-3af6-5f6dc915480f
# ╠═d16e91d0-a599-11ea-23f1-1963e0d078bb
# ╠═abd6ba19-16c4-4844-98bc-b1c7ee5716bd
# ╠═89a4e2a6-d7eb-11ea-3b22-fb24e94f2aa9
# ╠═a93ba57a-d7eb-11ea-26c8-1180ccb530d5
# ╠═4bd42778-d7ec-11ea-095b-8f80fca95014
# ╠═c2e29b42-7362-11eb-10e1-e1653cdd521e
