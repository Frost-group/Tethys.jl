include("Diagram.jl")
include("update.jl")
#include("measure.jl")
#include("zero_update.jl")
import Tethys
using Random
using Logging
using BenchmarkTools
using Profile

function insert_arc_benchmark(α, μ=-1.5, inserts=1, p=0)
    max_τ=30; max_order=500; mass=1; ω=1;
    diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)
    order=diagram.order
    m=mass
    m=mass
    α_squared=2pi*α*sqrt(2)
    for i in 1:inserts
        insert_arc!(diagram,order,mass, μ, ω, α_squared)
        order=diagram.order
        #swap_arc!(diagram)
        #scale!(diagram,order,m,μ,ω,1)
        
        #extend!(diagram)
        
    end
    #remove_arc2!(diagram,order,mass, μ, ω, α_squared)
end

begin
    t = @benchmark insert_arc_benchmark(1.0,-1.2,5000)
end

begin
    @profile insert_arc_benchmark(5.0,-6.0,100)
    Profile.print(format=:flat)
end

function full_benchmark(α, diagram, n_loop=1, p=0)
    max_τ=100.0; max_order=5000;
    n_hist=100000; sample_freq=2000;
    estimators = Estimators_Record(max_τ, max_order, sample_freq, n_loop, n_hist)
    diagram, estimators, variance = simulate!(diagram, estimators, true)
end

begin
    α=15.0
    loop_list = [1500]
    p=0.0
    max_τ=100.0; max_order=5000;
    n_hist=100000; sample_freq=2000;
    n_loop = loop_list[1]
    diagram = initialise_diagram(α, p, max_τ, max_order)
    estimators = Estimators_Record(max_τ, max_order, sample_freq, n_loop, n_hist)
    diagram, estimators, variance = simulate!(diagram, estimators, true)
end

begin
    BenchmarkTools.DEFAULT_PARAMETERS.seconds = 100
    t = @benchmark full_benchmark(15.0, diagram, 10)
end


function insert_arc_2_benchmark(α, μ, p=0)
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