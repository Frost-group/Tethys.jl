include("../src/new_measure.jl")
using Random
using LsqFit
using JLD2
using FFTW
using Logging
using Statistics
using Base.Threads

begin
    α=10.0; p=0.0
    max_τ=50.0; max_order=2000;
    n_loop=100; n_hist=100000; sample_freq=200;
    diagram = initialise_diagram(α, p, max_τ, max_order)
    estimators = Estimators_Record(max_τ, max_order, sample_freq, n_loop, n_hist)
    diagram, estimators, variance = simulate!(diagram, estimators, true)
end

begin
    estimators = Estimators_Record(max_τ, max_order, sample_freq, 100, n_hist)
    diagram, estimators, variance = simulate!(diagram, estimators, true)
end

begin
    α_list = collect(1.0:20.0)
    loop_list = [1000, 2000, 5000, 10000]
    p=0.0
    max_τ=30.0; max_order=2000;
    n_hist=100000; sample_freq=200;
    diagram_list = []
    estimators_list = []
    variance_list = []
    for α in α_list
        if α < 5.0
            n_loop = loop_list[1]
        elseif α < 10.0
            n_loop = loop_list[2]
        elseif α < 15.0
            n_loop = loop_list[3]
        else
            n_loop = loop_list[4]
        end
        diagram = initialise_diagram(α, p, max_τ, max_order)
        estimators = Estimators_Record(max_τ, max_order, sample_freq, n_loop, n_hist)
        diagram, estimators, variance = simulate!(diagram, estimators, false)
        push!(diagram_list, diagram)
        push!(estimators_list, estimators)
        push!(variance_list, variance)
    end
end

begin
    energy_list = []
    mass_list = []
    for estimator in estimators_list
        push!(energy_list, estimator.energy_mean[end])
        push!(mass_list, estimator.mass_mean[end])
    end
end