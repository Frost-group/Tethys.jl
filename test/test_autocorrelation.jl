include("../src/new_measure.jl")
using Random
using LsqFit
using JLD2
using FFTW
using Logging
using Statistics
using Base.Threads
using LaTeXStrings

begin
    α_list = collect(1.0:7.0)
    loop_list = [25]
    p=0.0
    max_τ=100.0; max_order=5000;
    n_hist=100000; sample_freq=1;
    diagram_listauto = []
    estimators_listauto = []
    auto_list = []
    #variance_listauto = []
    for α in α_list
        n_loop = loop_list[1]
        diagram = initialise_diagram(α, p, max_τ, max_order)
        estimators = Estimators_Record(max_τ, max_order, sample_freq, n_loop, n_hist)
        diagram, estimators, variance = simulate!(diagram, estimators, false)
        push!(diagram_listauto, diagram)
        push!(estimators_listauto, estimators)
        push!(auto_list, autocorrelation(estimators.energy_record,100000))
        #push!(variance_listmass, variance)
    end
end

begin
    α_list = collect(range(0.2, 5.0, length=25))
    loop_list = [10]
    p=0.0
    max_τ=100.0; max_order=5000;
    n_hist=100000; sample_freq=1;
    auto_list_mean = []
    auto_list_std = []
    #variance_listauto = []
    for α in α_list
        n_loop = loop_list[1]
        auto_vals = []
        for i in 1:50
            diagram = initialise_diagram(α, p, max_τ, max_order)
            estimators = Estimators_Record(max_τ, max_order, sample_freq, n_loop, n_hist)
            diagram, estimators, variance = simulate!(diagram, estimators, false)
            push!(auto_vals, autocorrelation(estimators.energy_record,100000))
        end
        push!(auto_list_mean, mean(auto_vals))
        push!(auto_list_std, std(auto_vals))
        #push!(variance_listmass, variance)
    end
end

begin
    plot_font = "Computer Modern"
    default(fontfamily=plot_font, dpi = 500)
    scatter(collect(range(0.2, 5.0, length=25)),auto_list_mean,yerr=auto_list_std./sqrt(50), ylab = L"$\langle R_{t} \rangle$",
    xlab = L"\alpha", framestyle = :box, xtickfontsize=10,ytickfontsize=10, gridlinewidth = 1.0,
    guidefontsize = 14, linewidth=1.5, markerstrokewidth=1.0, legendfontsize = 10, legend = false,
    markersize=2.5, markeralpha=1)

    #savefig("C:/Users/60125/OneDrive - Imperial College London/Y4 Masters Project/Graphics/autocorrelation")
end

begin
    α_list = collect(range(0.2, 5.0, length=25))
    loop_list = [25]
    p=0.0
    max_τ=100.0; max_order=5000;
    n_hist=100000; sample_freq=200;
    order_mean = []
    order_std = []
    #variance_listauto = []
    for α in α_list
        n_loop = loop_list[1]
        diagram = initialise_diagram(α, p, max_τ, max_order)
        estimators = Estimators_Record(max_τ, max_order, sample_freq, n_loop, n_hist)
        diagram, estimators, variance = simulate!(diagram, estimators, false)
        push!(order_mean, argmax(estimators.order_box))
        push!(order_std, std(estimators.order_box))
        #push!(variance_listmass, variance)
    end
end

begin
    plot_font = "Computer Modern"
    default(fontfamily=plot_font, dpi = 500)
    scatter(collect(range(0.2, 5.0, length=25)),order_mean, yerr=order_std, ylab = L"$\langle N \rangle$",
    xlab = L"\alpha", framestyle = :box, xtickfontsize=10,ytickfontsize=10, gridlinewidth = 1.0,
    guidefontsize = 14, linewidth=1.5, markerstrokewidth=1.0, legendfontsize = 10, legend = false,
    markersize=2.5, markeralpha=1)

    #savefig("C:/Users/60125/OneDrive - Imperial College London/Y4 Masters Project/Graphics/order")
end

begin
    α=11.0
    loop_list = [1000]
    p=0.0
    max_τ=100.0; max_order=5000;
    n_hist=100000; sample_freq=2000;
    n_loop = loop_list[1]
    diagram = initialise_diagram(α, p, max_τ, max_order)
    estimators = Estimators_Record(max_τ, max_order, sample_freq, n_loop, n_hist)
    diagram, estimators, variance = simulate!(diagram, estimators, false)

    ar = estimators.mass_mean
    plot(ar)
end

begin
    α=11.0
    loop_list = [1000]
    p=0.0
    max_τ=100.0; max_order=5000;
    n_hist=100000; sample_freq=2000;
    n_loop = loop_list[1]
    diagram = initialise_diagram(α, p, max_τ, max_order)
    estimators = Estimators_Record(max_τ, max_order, sample_freq, n_loop, n_hist)
    diagram, estimators, variance = simulate!(diagram, estimators, true)

    br = estimators.mass_mean
    plot!(br)
end

begin
    plot_font = "Computer Modern"
    default(fontfamily=plot_font, dpi = 500)
    plot(ar[10:end], ylab = L"$\langle R_{insert} \rangle$", label=L"\mathrm{w/o}\:\:swap", color=:red,
    xlab = L"Samples", framestyle = :box, xtickfontsize=10,ytickfontsize=10, gridlinewidth = 0.5,
    guidefontsize = 12, linewidth=1.2, markerstrokewidth=1.0, legendfontsize = 11, foreground_color_legend = nothing,
    markersize=2.5, markeralpha=1, extra_kwargs=Dict(:subplot=>Dict(:legend_hfactor=>1.25)))
    plot!(br[10:end], ylab = L"$\langle R_{insert} \rangle$", label=L"\mathrm{with}\:\:swap", color=:green,
    xlab = "Samples", framestyle = :box, xtickfontsize=10,ytickfontsize=10, gridlinewidth = 0.5,
    guidefontsize = 12, linewidth=1.2, markerstrokewidth=1.0, legendfontsize = 11, foreground_color_legend = nothing,
    markersize=2.5, markeralpha=1, extra_kwargs=Dict(:subplot=>Dict(:legend_hfactor=>1.25)))

    #savefig("C:/Users/60125/OneDrive - Imperial College London/Y4 Masters Project/Graphics/acceptance_ratio")
end