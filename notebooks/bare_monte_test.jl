### A Pluto.jl notebook ###
# v0.9.2

using Markdown

# ╔═╡ 6d6a567a-a593-11ea-2070-53d683c984d6
# boot local version, under Revise
begin
	using Revise
	using Pkg
	Pkg.activate("/Users/jarvist/REPOS/Tethys.jl/") # ugly full path...
	using Tethys
	
	using Plots
	using StatsPlots
end

# ╔═╡ 9a8c77a0-a593-11ea-09e0-5177021a936e
begin
	MAX_ORDER=2
	diag=Tethys.Diag(0,1.0,MAX_ORDER,[0.1,0.2,1.0 , 0.3,0.4,1.0], [0,0,0,0,0,0])

	GF1=Tethys.BareExpansion(diag)
end

# ╔═╡ ddcd5f98-a593-11ea-1b0f-6f6067269fc4
GF2=Tethys.Monte!(diag)

# ╔═╡ ed0ad5da-a593-11ea-0d19-1fa974170c9d
begin
	GF=Tethys.BareExpansion(diag)
	MC=[]
	for i in 1:100
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

# ╔═╡ Cell order:
# ╠═6d6a567a-a593-11ea-2070-53d683c984d6
# ╠═9a8c77a0-a593-11ea-09e0-5177021a936e
# ╠═ddcd5f98-a593-11ea-1b0f-6f6067269fc4
# ╠═ed0ad5da-a593-11ea-0d19-1fa974170c9d
# ╠═6bef50b2-a594-11ea-3ab0-dd19b06da250
# ╠═6f29880e-a594-11ea-3a1e-b18693818f43
# ╠═e8d22018-a593-11ea-2d09-c921536e9287
