include("Diagram.jl")
include("update.jl")
include("measure.jl")

using Profile
using Logging
using BenchmarkTools

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
    μ_list = -α.-1.26*(α./10).^2 .+[-3.0:0.4:3.0;]
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
    directory = "D://data"
    diag_order = []
    α = 1
    μ_list = -α.-1.26*(α./10).^2 .+[-3.0:0.4:3.0;]
    p=0; max_τ=30; max_order=500; mass=1; ω=1;
    for μ in μ_list
        @info "Loading μ = "*string(μ)
        n_loop = 2000
        store_data = false

        hist=Hist_Record(300,max_τ,max_order)
        diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)

        diagram,hist,green_record,zero_record,normalized_data,bin_variance  = hist_measure!(diagram,hist,directory,n_loop,store_data);
        append!(diag_order, diagram.order)

        @info "Diagram order" diagram.order
    end
end

begin
    plot(μ_list,diag_order,xlabel="μ",ylabel="Diagram order",label="α = "*string(α))
end

begin
    directory = "F://data"
    directory_list = readdir(directory; join=true)
    for dir in directory_list
        α = parse(Float64, split(dir, "=")[2])
        @info "Loading α = "*string(α)
        timings = []
        μ_list = -α.-1.26*(α./10).^2 .+[-3.0:0.4:3.0;]
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

begin
    n_loop=200
    num_mea=1; regime=Diff_more();
    p=0; max_τ=30; max_order=500; mass=1; ω=1;


    linear(t, p) = p[1].-p[2].*t
    bin_width=max_τ/300
    min_time=Int(div(5,bin_width,RoundUp))
    max_time=Int(div(12,bin_width,RoundUp))

    α=0.1
    μ=1
    directory = "F://data"
    store_data = false
    hist=Hist_Record(300,max_τ,max_order)
    diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)
    diagram,hist,green_record,zero_record,green_func,variance=hist_measure!(diagram,hist,directory,n_loop,store_data)

    time_points=hist.time_points[min_time:max_time]

    statis=sum(green_func[i,:] for i in 1:max_order+1)
    y=log.(statis)[min_time:max_time]
    w=1 ./(variance[min_time:max_time]./(statis[min_time:max_time]).^2)

    p0=[0,(-α-1.26*(α/10)^2-μ)]
    fit = curve_fit(linear, time_points, y, w, p0)
    errors=standard_errors(fit)


    z0=exp(fit.param[1])
    plot(time_points,y)
    display(plot!(time_points,linear(time_points,fit.param),
    xlabel="τ",ylabel="log(green)",title="α="*string(α)*",k="*string(0)*",μ="*string(μ)))
    println("energy:",fit.param[2]+μ)
    println("perturb:",-α-1.26*(α/10)^2)
end


