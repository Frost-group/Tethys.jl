include("Diagram.jl")
include("update.jl")
include("measure.jl")

using Profile
using Logging

begin
    p=0; max_τ=30; max_order=500; mass=1; ω=1;
    α = 1.1
    μ = -1.2
    n_loop = 10
    store_data = false

    hist=Hist_Record(300,max_τ,max_order)
    diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)

    @profile hist_measure!(diagram,hist,"D://data",n_loop,store_data);
    Profile.print(format=:flat)

end

begin
    p=0; max_τ=30; max_order=500; mass=1; ω=1;
    α = 1.1
    μ = -1.2
    n_loop = 10
    store_data = false

    hist=Hist_Record(300,max_τ,max_order)
    diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)

    stats = @timed hist_measure!(diagram,hist,"D://data",n_loop,store_data);

    @info "Time per sweep w/o data collection" stats.time/n_loop
    @info "Memory allocated" stats.bytes
end

begin
    p=0; max_τ=30; max_order=500; mass=1; ω=1;
    α = 1.1
    μ = -1.2
    n_loop = 10
    store_data = true

    hist=Hist_Record(300,max_τ,max_order)
    diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)

    stats = @timed hist_measure!(diagram,hist,"D://data",n_sweep,store_data);

    @info "Time per sweep with data collection" stats.time/n_loop
    @info "Memory allocated" stats.bytes
end

begin
    p=0; max_τ=30; max_order=500; mass=1; ω=1;
    α = 1.1
    μ = -1.2
    n_loop = 10
    store_data = false
    n_hist = 10000

    hist=Hist_Record(300,max_τ,max_order)
    diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)

    stats = @timed hist_measure!(diagram,hist,"D://data",n_loop,store_data,n_hist);

    @info "Time per sweep with 10000 histograms sampled" stats.time/n_loop
    @info "Memory allocated" stats.bytes
end

begin
    p=0; max_τ=30; max_order=500; mass=1; ω=1;
    α = 1.1
    μ = -1.2
    n_loop = 10
    store_data = false
    n_hist = 1000000

    hist=Hist_Record(300,max_τ,max_order)
    diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)

    stats = @timed hist_measure!(diagram,hist,"D://data",n_loop,store_data,n_hist);

    @info "Time per sweep with 1000000 histograms sampled" stats.time/n_loop
    @info "Memory allocated" stats.bytes
end

begin
    p=0; max_τ=30; max_order=500; mass=1; ω=1;
    α = 1.1
    μ = -1.2
    n_loop = 10
    store_data = false
    n_hist = 10000

    hist=Hist_Record(300,max_τ,max_order)
    diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)

    diagram,hist,green_record,zero_record,normalized_data,bin_variance = hist_measure!(diagram,hist,"D://data",n_loop,store_data,n_hist);
end

begin
    directory = "D://data"
    timings = []
    μ_list = α.-1.26*(α./10).^2 .+[-1.0:0.1:1.0;]
    p=0; max_τ=30; max_order=500; mass=1; ω=1;
    α = 0.1
    for μ in μ_list
        n_loop = 100
        store_data = false

        hist=Hist_Record(300,max_τ,max_order)
        diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)

        stats = @timed hist_measure!(diagram,hist,directory,n_loop,store_data);
        append!(timings, stats.time/n_loop)

        @info "Time per sweep w/o data collection" stats.time/n_loop
        @info "Memory allocated" stats.bytes
    end
end

begin
    plot(μ_list,timings,xlabel="μ",ylabel="Seconds per sweep",label="α = "*string(α))
end

begin
    directory = "F://data"
    directory_list = readdir(directory; join=true)
    for dir in directory_list
        α = parse(Float64, split(dir, "=")[2])
        @info "Loading α = "*string(α)
        timings = []
        μ_list = -α.-1.26*(α./10).^2 .+[-1.0:0.1:1.0;]
        p=0; max_τ=30; max_order=500; mass=1; ω=1;
        for μ in μ_list
            @info "Loading μ = "*string(μ)
            n_loop = 100
            store_data = false

            hist=Hist_Record(300,max_τ,max_order)
            diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)

            stats = @timed hist_measure!(diagram,hist,directory,n_loop,store_data);
            append!(timings, stats.time/n_loop)

            @info "Time per sweep w/o data collection" stats.time/n_loop
            @info "Memory allocated" stats.bytes
        end
        display(plot(μ_list,timings,xlabel="μ",ylabel="Seconds per sweep",label="α = "*string(α)))
    end
end
