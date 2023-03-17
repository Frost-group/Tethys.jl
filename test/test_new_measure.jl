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
    α=1.0; p=1.5
    max_τ=100.0; max_order=5000;
    n_loop=100; n_hist=100000; sample_freq=200;
    diagram = initialise_diagram(α, p, max_τ, max_order)
    estimators = Estimators_Record(max_τ, max_order, sample_freq, n_loop, n_hist)
    diagram, estimators, variance = simulate!(diagram, estimators, false)
end

begin
    estimators = Estimators_Record(max_τ, max_order, sample_freq, 100, n_hist)
    diagram, estimators, variance = simulate!(diagram, estimators, true)
end

begin
    α_list = collect(1.0:10.0)
    loop_list = [1000, 2000, 5000, 10000]
    p=0.0
    max_τ=50.0; max_order=2000;
    n_hist=100000; sample_freq=200;
    diagram_list = []
    estimators_list = []
    variance_list = []
    for α in α_list
        if α < 5.0
            n_loop = loop_list[1]
        elseif α < 11.0
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

begin
    α_list = collect(11.0:16.0)
    loop_list = [2000]
    p=0.0
    max_τ=100.0; max_order=5000;
    n_hist=100000; sample_freq=200;
    diagram_list2 = []
    estimators_list2 = []
    variance_list2 = []
    for α in α_list
        n_loop = loop_list[1]
        diagram = initialise_diagram(α, p, max_τ, max_order)
        estimators = Estimators_Record(max_τ, max_order, sample_freq, n_loop, n_hist)
        diagram, estimators, variance = simulate!(diagram, estimators, true)
        push!(diagram_list2, diagram)
        push!(estimators_list2, estimators)
        push!(variance_list2, variance)
    end
end

begin
    energy_list2 = []
    mass_list2 = []
    for estimator in estimators_list2
        push!(energy_list2, estimator.energy_mean[end])
        push!(mass_list2, estimator.mass_mean[end])
    end
end

begin
    α=1.0
    max_τ=100.0; max_order=5000;
    n_loop=500; n_hist=100000; sample_freq=200;
    k_list = collect(range(0, 2.0, length=25))
    diagram_klist = []
    estimators_klist = []
    variance_klist = []
    for k in k_list
        diagram = initialise_diagram(α, k, max_τ, max_order)
        estimators = Estimators_Record(max_τ, max_order, sample_freq, n_loop, n_hist)
        diagram, estimators, variance = simulate!(diagram, estimators, false)
        push!(diagram_klist, diagram)
        push!(estimators_klist, estimators)
        push!(variance_klist, variance)
    end
end

begin
    α=1.0
    max_τ=100.0; max_order=5000;
    n_loop=1000; n_hist=100000; sample_freq=200;
    k_list = collect(range(1.5, 1.9, length=10))
    diagram_klist2 = []
    estimators_klist2 = []
    variance_klist2 = []
    for k in k_list
        diagram = initialise_diagram(α, k, max_τ, max_order)
        estimators = Estimators_Record(max_τ, max_order, sample_freq, n_loop, n_hist)
        diagram, estimators, variance = simulate!(diagram, estimators, false)
        push!(diagram_klist2, diagram)
        push!(estimators_klist2, estimators)
        push!(variance_klist2, variance)
    end
end

begin
    energy_klist2 = []
    mass_klist2 = []
    variance_klist2 = []
    for estimator in estimators_klist2
        push!(energy_klist2, estimator.energy_mean[end])
        push!(mass_klist2, estimator.mass_mean[end])
        push!(variance_klist2, std(estimator.energy_mean[100000:end])^2)
    end
end

begin
    α=2.0
    max_τ=100.0; max_order=5000;
    n_loop=500; n_hist=100000; sample_freq=200;
    k_list = collect(range(0, 2.0, length=25))
    diagram_k2list = []
    estimators_k2list = []
    variance_k2list = []
    for k in k_list
        diagram = initialise_diagram(α, k, max_τ, max_order)
        estimators = Estimators_Record(max_τ, max_order, sample_freq, n_loop, n_hist)
        diagram, estimators, variance = simulate!(diagram, estimators, false)
        push!(diagram_k2list, diagram)
        push!(estimators_k2list, estimators)
        push!(variance_k2list, variance)
    end
end

begin
    energy_klist = []
    mass_klist = []
    variance_klist = []
    for estimator in estimators_klist
        push!(energy_klist, estimator.energy_mean[end])
        push!(mass_klist, estimator.mass_mean[end])
        push!(variance_klist, std(estimator.energy_mean[100000:end])^2)
    end
end


begin
    quadratic(t, p) = p[1].+p[2].*t.+p[3].*t.*t
    p0=[0.0,-1.0,1.0]
    w=1 ./sqrt.(variance_klist[1:10])
    fit = curve_fit(quadratic, k_list[1:10], energy_klist[1:10], w,p0)
    plot_font = "Computer Modern"
    default(fontfamily=plot_font, dpi = 500)

    plot(k_list,  energy_klist, yerr=sqrt.(variance_klist), label = L"$\alpha=1$", ylab = L"$E_{0}(k)$",
    xlab = L"k", framestyle = :box, xtickfontsize=10,ytickfontsize=10, gridlinewidth = 0.0,
    guidefontsize = 14, linewidth=1.5, markerstrokewidth=1.0, legendfontsize = 10, legend=:topleft,
    foreground_color_legend = nothing, extra_kwargs=Dict(:subplot=>Dict(:legend_hfactor=>1.25)))
    plot!(k_list, quadratic(k_list, fit.param), label = "Fit", linestyle=:dash)
    hline!([energy_klist[1]+1], linestyle=:dot, label = "Continuum Edge")
    #plot!(alpha_plot_list, VMC_plot_list, label = "Feynman")


    #savefig("C:/Users/60125/OneDrive - Imperial College London/Y4 Masters Project/Graphics/dispersion_alpha_1")
end

begin
    energy_klist_total = vcat(energy_klist[1:18], energy_klist2[1:8])
    k_list_total = vcat(collect(range(0, 2.0, length=25))[1:18], collect(range(1.5, 1.9, length=10))[1:8])
    variance_klist_total = vcat(variance_klist[1:18], variance_klist2[1:8])
    quadratic(t, p) = p[1].+p[2].*t.+p[3].*t.*t
    p0=[0.0,-1.0,1.0]
    w=1 ./sqrt.(variance_klist_total[1:8])
    fit = curve_fit(quadratic, k_list_total[1:8], energy_klist_total[1:8], w,p0)
    plot_font = "Computer Modern"
    default(fontfamily=plot_font, dpi = 500)

    plot(k_list_total,  energy_klist_total, yerr=sqrt.(variance_klist_total), label = L"$\alpha=1$", ylab = L"$E_{0}(k)$",
    xlab = L"k", framestyle = :box, xtickfontsize=10,ytickfontsize=10, gridlinewidth = 0.0,
    guidefontsize = 14, linewidth=2.0, markerstrokewidth=1.0, legendfontsize = 10, legend=:topleft,
    foreground_color_legend = nothing, extra_kwargs=Dict(:subplot=>Dict(:legend_hfactor=>1.25)))
    plot!(k_list_total, quadratic(k_list_total, fit.param), label = "Parabolic fit", linestyle=:dash, linecolor=:red)
    hline!([energy_klist_total[1]+1], linestyle=:dot, label = "Continuum edge", linecolor=:grey)
    #plot!(alpha_plot_list, VMC_plot_list, label = "Feynman")


    #savefig("C:/Users/60125/OneDrive - Imperial College London/Y4 Masters Project/Graphics/dispersion_alpha_1")
end

begin
    energy_k2list = []
    mass_k2list = []
    variance_k2list = []
    for estimator in estimators_k2list
        push!(energy_k2list, estimator.energy_mean[end])
        push!(mass_k2list, estimator.mass_mean[end])
        push!(variance_k2list, std(estimator.energy_mean[100000:end])^2)
    end
end

begin
    α=2.0
    max_τ=100.0; max_order=5000;
    n_loop=1000; n_hist=100000; sample_freq=200;
    k_list = collect(range(2.0, 2.5, length=10))
    diagram_k2list2 = []
    estimators_k2list2 = []
    variance_k2list2 = []
    for k in k_list
        diagram = initialise_diagram(α, k, max_τ, max_order)
        estimators = Estimators_Record(max_τ, max_order, sample_freq, n_loop, n_hist)
        diagram, estimators, variance = simulate!(diagram, estimators, false)
        push!(diagram_k2list2, diagram)
        push!(estimators_k2list2, estimators)
        push!(variance_k2list2, variance)
    end
end

begin
    energy_k2list2 = []
    mass_k2list2 = []
    variance_k2list2 = []
    for estimator in estimators_k2list2
        push!(energy_k2list2, estimator.energy_mean[end])
        push!(mass_k2list2, estimator.mass_mean[end])
        push!(variance_k2list2, std(estimator.energy_mean[100000:end])^2)
    end
end

begin
    energy_klist_total = vcat(energy_k2list, energy_k2list2[2:8])
    k_list_total = vcat(collect(range(0, 2.0, length=25)), collect(range(2.0, 2.5, length=10))[2:8])
    variance_klist_total = vcat(variance_k2list, variance_k2list2[2:8])
    quadratic(t, p) = p[1].+p[2].*t.+p[3].*t.*t
    p0=[0.0,-1.0,1.0]
    w=1 ./sqrt.(variance_klist_total[1:8])
    fit = curve_fit(quadratic, k_list_total[1:8], energy_klist_total[1:8], w,p0)
    plot_font = "Computer Modern"
    default(fontfamily=plot_font, dpi = 500)

    plot(k_list_total,  energy_klist_total, yerr=sqrt.(variance_klist_total), label = L"$\alpha=2$", ylab = L"$E_{0}(k)$",
    xlab = L"k", framestyle = :box, xtickfontsize=10,ytickfontsize=10, gridlinewidth = 0.0,
    guidefontsize = 14, linewidth=2.0, markerstrokewidth=1.0, legendfontsize = 10, legend=:topleft,
    foreground_color_legend = nothing, extra_kwargs=Dict(:subplot=>Dict(:legend_hfactor=>1.25)))
    plot!(k_list_total, quadratic(k_list_total, fit.param), label = "Parabolic fit", linestyle=:dash, linecolor=:red)
    hline!([energy_klist_total[1]+1], linestyle=:dot, label = "Continuum edge", linecolor=:grey)
    #plot!(alpha_plot_list, VMC_plot_list, label = "Feynman")


    #savefig("C:/Users/60125/OneDrive - Imperial College London/Y4 Masters Project/Graphics/dispersion_alpha_2")
end
