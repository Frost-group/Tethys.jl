include("Diagram.jl")
include("update.jl")
include("measure.jl")

using Logging
using BenchmarkTools

function insert_arc_benchmark(α, μ, inserts=1, p=0)
    max_τ=30; max_order=500; mass=1.0; ω=1.0;
    diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)
    order=diagram.order
    m=mass
    p=diagram.p
    dispersion_val = norm(p)^2/(2m)-μ
    α_squared=2pi*α*sqrt(2)
    for i in 1:inserts
        insert_arc!(diagram,order,mass, μ, ω, α_squared)
        extend!(diagram, dispersion_val)
        order = diagram.order
    end
end

begin
    t = @benchmark insert_arc_benchmark(5.0, -6.0,1)
end

function extend_benchmark(α, μ, p=0)
    max_τ=30; max_order=500; mass=1; ω=1;
    regime=Diff_more()
    diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)
    extend!(diagram)
end

begin
    t = @benchmark extend_benchmark(5.0, -6)
end


begin
    n_loop=10
    num_mea=1; regime=Diff_more();
    p=0; max_τ=30; max_order=500; mass=1; ω=1;


    linear(t, p) = p[1].-p[2].*t
    bin_width=max_τ/300
    min_time=Int(div(5,bin_width,RoundUp))
    max_time=Int(div(12,bin_width,RoundUp))

    α=5
    μ=-6
    directory = "F://data"
    store_data = false
    hist=Hist_Record(300,max_τ,max_order)
    diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)
    diagram,hist,green_record,zero_record,green_func,variance=hist_measure!(diagram,hist,directory,n_loop,store_data)
    
    @info "Diagram Order: "*string(diagram.order)
    t = @benchmark insert_arc!(diagram,regime)
end

function measure_loop(α, μ, p=0)
    n_loop=100
    num_mea=1; regime=Diff_more();
    max_τ=30; max_order=500; mass=1; ω=1;

    directory = "F://data"
    store_data = false
    hist=Hist_Record(300,max_τ,max_order)
    diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)
    diagram,hist,green_record,zero_record,green_func=hist_measure!(diagram,hist,directory,n_loop,store_data,1000)
    
end

begin
    Random.seed!(1234)
    BenchmarkTools.DEFAULT_PARAMETERS.seconds = 5
    t = @benchmark measure_loop(1, -1.2)
end

begin
    histogram([t2.times])
end