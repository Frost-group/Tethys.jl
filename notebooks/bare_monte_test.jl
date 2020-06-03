### A Pluto.jl notebook ###
# v0.9.2

using Markdown

# ╔═╡ b6e035fc-a5a9-11ea-04b4-f9c64fe7219d
begin
	using Plots
	using StatsPlots
end

# ╔═╡ 6d6a567a-a593-11ea-2070-53d683c984d6
# boot local version, under Revise
begin
	using Revise
	using Pkg
	Pkg.activate("/Users/jarvist/REPOS/Tethys.jl/") # ugly full path...
	using Tethys
	
end

# ╔═╡ 9a8c77a0-a593-11ea-09e0-5177021a936e
begin
	MAX_ORDER=2
	diag1=Tethys.Diag(0,1.0,MAX_ORDER,[0.1,0.2,1.0 , 0.3,0.4,1.0], [0.0,0.0,0.0,0.0,0.0,0.0])
	GF1=Tethys.BareExpansion(diag1)
end

# ╔═╡ fa04d3d6-a5ab-11ea-1be2-fb7c20cdeb52


# ╔═╡ 32209058-a5ab-11ea-0306-5d121ab3c45b


# ╔═╡ 2399f97a-a5ab-11ea-11cc-0934ede99f55


# ╔═╡ 6c3748a8-a5a7-11ea-2490-d5e6dce272a8
diag1

# ╔═╡ ddcd5f98-a593-11ea-1b0f-6f6067269fc4
GF2=Tethys.Monte!(diag1)

# ╔═╡ 69b81526-a5a7-11ea-272f-21121dbbfbb9
diag1

# ╔═╡ ed0ad5da-a593-11ea-0d19-1fa974170c9d
begin
	diag=Tethys.Diag(0,1.0,MAX_ORDER,[0.1,0.2,1.0 , 0.3,0.4,1.0], [0,0,0,0,0,0])
	MC=[]
	for i in 1:100_000
		GF=Tethys.Monte!(diag)
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
Tethys.BareExpansion(diag)

# ╔═╡ b9e77af0-a5a7-11ea-3af8-5bd668015649
if (1==2)==false
	a=2
end

# ╔═╡ e74ffc6a-a5a7-11ea-02a9-7373c4c62c17


# ╔═╡ c61425d0-a5a7-11ea-0268-f51f1440167e


# ╔═╡ c2d4a91c-a5a7-11ea-2630-db855e65d7bf


# ╔═╡ c0859b30-a5a7-11ea-0634-c98eb52c37e8


# ╔═╡ Cell order:
# ╠═b6e035fc-a5a9-11ea-04b4-f9c64fe7219d
# ╠═6d6a567a-a593-11ea-2070-53d683c984d6
# ╠═9a8c77a0-a593-11ea-09e0-5177021a936e
# ╠═fa04d3d6-a5ab-11ea-1be2-fb7c20cdeb52
# ╠═32209058-a5ab-11ea-0306-5d121ab3c45b
# ╠═2399f97a-a5ab-11ea-11cc-0934ede99f55
# ╠═6c3748a8-a5a7-11ea-2490-d5e6dce272a8
# ╠═ddcd5f98-a593-11ea-1b0f-6f6067269fc4
# ╠═69b81526-a5a7-11ea-272f-21121dbbfbb9
# ╠═ed0ad5da-a593-11ea-0d19-1fa974170c9d
# ╠═6bef50b2-a594-11ea-3ab0-dd19b06da250
# ╠═6f29880e-a594-11ea-3a1e-b18693818f43
# ╠═e8d22018-a593-11ea-2d09-c921536e9287
# ╠═d08f9624-a599-11ea-3af6-5f6dc915480f
# ╠═d16e91d0-a599-11ea-23f1-1963e0d078bb
# ╠═b9e77af0-a5a7-11ea-3af8-5bd668015649
# ╠═e74ffc6a-a5a7-11ea-02a9-7373c4c62c17
# ╠═c61425d0-a5a7-11ea-0268-f51f1440167e
# ╠═c2d4a91c-a5a7-11ea-2630-db855e65d7bf
# ╠═c0859b30-a5a7-11ea-0634-c98eb52c37e8
