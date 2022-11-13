include("Diagram.jl")
include("update.jl")
include("measure.jl")

using Logging
using BenchmarkTools

function insert_arc_benchmark(α, μ, p=0)
    max_τ=30; max_order=500; mass=1; ω=1;
    regime=Diff_more()
    diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)
    insert_arc!(diagram,regime)
end

begin
    t = @benchmark insert_arc_benchmark(5.0, -6)
end

function testing_arr(arr)
    q,w,e,r = arr

end
begin
    t = @benchmark testing_arr([1,2,3,4])
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
    n_loop=10
    num_mea=1; regime=Diff_more();
    max_τ=30; max_order=500; mass=1; ω=1;


    linear(t, p) = p[1].-p[2].*t
    bin_width=max_τ/300
    min_time=Int(div(5,bin_width,RoundUp))
    max_time=Int(div(12,bin_width,RoundUp))

    directory = "F://data"
    store_data = false
    hist=Hist_Record(300,max_τ,max_order)
    diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)
    diagram,hist,green_record,zero_record,green_func,variance=hist_measure!(diagram,hist,directory,n_loop,store_data,1000)
    
end

begin
    t = @benchmark measure_loop(5.0, -6)
end
