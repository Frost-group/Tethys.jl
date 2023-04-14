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
    α=7.0; p=0.0
    max_τ=1000.0; max_order=5000;
    n_loop=500; n_hist=100000; sample_freq=200;
    diagram = initialise_diagram(α, p, max_τ, max_order)
    estimators = Estimators_Record(max_τ, max_order, sample_freq, n_loop, n_hist)
    diagram, estimators, variance = simulate!(diagram, estimators,false)
end

begin
    plot_font = "Computer Modern"
    default(fontfamily=plot_font, dpi = 500)

    plot(estimators.energy_mean[10:end], label = "DiagMC", ylab = L"$\langle E_{0}(\alpha=1)\rangle$",
    xlab = "Samples", framestyle = :box, xtickfontsize=8,ytickfontsize=9, color=:red,
    xguidefontsize = 12, linewidth=1.5, markerstrokewidth=1.0, legendfontsize = 11, legend=:topright,
    foreground_color_legend = nothing, extra_kwargs=Dict(:subplot=>Dict(:legend_hfactor=>1.25)),
    yguidefontsize=14)
    hline!([VMC_energy(1)], linestyle=:dash, label = "Feynman", linecolor=:black, linewidth=1.0)
    #savefig("C:/Users/60125/OneDrive - Imperial College London/Y4 Masters Project/Graphics/mean_energy_evolution_1")
end

begin
    plot_font = "Computer Modern"
    default(fontfamily=plot_font, dpi = 500)

    plot!(estimators.order_box[1:200], label = L"\alpha=7", ylab = "Order Distribution",
    xlab = L"N_{arcs}", framestyle = :box, xtickfontsize=8,ytickfontsize=9,
    xguidefontsize = 12, linewidth=1.0, markerstrokewidth=1.0, legendfontsize = 11, legend=:topright,
    foreground_color_legend = nothing, extra_kwargs=Dict(:subplot=>Dict(:legend_hfactor=>1.25)),
    yguidefontsize=10)
    #savefig("C:/Users/60125/OneDrive - Imperial College London/Y4 Masters Project/Graphics/diagram_order")
end


begin
    histogram(estimators.mass_mean[10000:end], label = "No Swap", xlab = "Mean line length", ylab = "Frequency")
end

begin
    histogram!(estimators.mass_mean[10000:end], label = "Swap", xlab = "Mean line length", ylab = "Frequency", alpha = 0.3)
    #savefig("C:/Users/60125/OneDrive - Imperial College London/Y4 Masters Project/Graphics/mean_line_length_15")
end

begin
    diagram = set_α!(diagram, 15.0)
    estimators = Estimators_Record(max_τ, max_order, sample_freq, 200, n_hist)
    diagram, estimators, variance = simulate!(diagram, estimators, true)
end


begin
    norm_list = []
    for i in 1:20000
        phi = rand(Uniform(0,pi*2))
        costheta = rand(Uniform(-1,1))
        theta = acos(costheta)
        x = sin(theta)*cos(phi)
        y = sin(theta)*sin(phi)
        z = cos(theta)
        q = abs(rand(Normal(0,sqrt(1/1)))).*[x,y,z]
        push!(norm_list,norm(q)^2)
    end
end

begin
    ave = [0.0, 0.0, 0.0]
    cov = [1 0 0; 0 1 0; 0 0 1]
    d = MvNormal(ave, cov)
    norm_list = []
    for i in 1:20000
        x = rand(d)
        #x = rand(Normal(0,1),3)
        push!(norm_list,norm(x)^2)
    end
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
    variance_list = []
    for estimator in estimators_list
        push!(energy_list, estimator.energy_mean[end])
        push!(mass_list, estimator.mass_mean[end])
        push!(variance_list, std(estimator.energy_mean[100000:end])^2)
    end
end

#swap
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
    variance_list2 = []
    for estimator in estimators_list2
        push!(energy_list2, estimator.energy_mean[end])
        push!(mass_list2, estimator.mass_mean[end])
        push!(variance_list2, std(estimator.energy_mean[100000:end])^2)
    end
end

#resample
begin
    α_list = collect(11.0:16.0)
    loop_list = [2000]
    p=0.0
    max_τ=100.0; max_order=5000;
    n_hist=100000; sample_freq=200;
    diagram_list3 = []
    estimators_list3 = []
    variance_list3 = []
    for α in α_list
        n_loop = loop_list[1]
        diagram = initialise_diagram(α, p, max_τ, max_order)
        estimators = Estimators_Record(max_τ, max_order, sample_freq, n_loop, n_hist)
        diagram, estimators, variance = simulate!(diagram, estimators, false)
        push!(diagram_list3, diagram)
        push!(estimators_list3, estimators)
        push!(variance_list3, variance)
    end
end

begin
    energy_list3 = []
    mass_list3 = []
    variance_list3 = []
    for estimator in estimators_list3
        push!(energy_list3, estimator.energy_mean[end])
        push!(mass_list3, estimator.mass_mean[end])
        push!(variance_list3, std(estimator.energy_mean[100000:end])^2)
    end
end

#swap
begin
    α_list = collect(16.0:20.0)
    loop_list = [4000]
    p=0.0
    max_τ=100.0; max_order=5000;
    n_hist=100000; sample_freq=200;
    diagram_list4 = []
    estimators_list4 = []
    variance_list4 = []
    for α in α_list
        n_loop = loop_list[1]
        diagram = initialise_diagram(α, p, max_τ, max_order)
        estimators = Estimators_Record(max_τ, max_order, sample_freq, n_loop, n_hist)
        diagram, estimators, variance = simulate!(diagram, estimators, true)
        push!(diagram_list4, diagram)
        push!(estimators_list4, estimators)
        push!(variance_list4, variance)
    end
end

begin
    energy_list4 = []
    mass_list4 = []
    variance_list4 = []
    for estimator in estimators_list4
        push!(energy_list4, estimator.energy_mean[end])
        push!(mass_list4, estimator.mass_mean[end])
        push!(variance_list4, std(estimator.energy_mean[100000:end])^2)
    end
end


begin
    α_list = collect(21.0:25.0)
    loop_list = [6000]
    p=0.0
    max_τ=100.0; max_order=5000;
    n_hist=100000; sample_freq=200;
    diagram_list5 = []
    estimators_list5 = []
    variance_list5 = []
    for α in α_list
        n_loop = loop_list[1]
        diagram = initialise_diagram(α, p, max_τ, max_order)
        estimators = Estimators_Record(max_τ, max_order, sample_freq, n_loop, n_hist)
        diagram, estimators, variance = simulate!(diagram, estimators, true)
        push!(diagram_list5, diagram)
        push!(estimators_list5, estimators)
        push!(variance_list5, variance)
    end
end

begin
    energy_list5 = []
    mass_list5 = []
    variance_list5 = []
    for estimator in estimators_list5
        push!(energy_list5, estimator.energy_mean[end])
        push!(mass_list5, estimator.mass_mean[end])
        push!(variance_list5, std(estimator.energy_mean[500000:end])^2)
    end
end


begin
    energy_list_total = vcat(energy_list, energy_list3)
    variance_list_total = vcat(variance_list, variance_list3)
    plot_font = "Computer Modern"
    default(fontfamily=plot_font, dpi = 500)

    plot(collect(1.0:16.0),  energy_list_total, yerr=sqrt.(variance_list_total), label = "DiagMC", ylab = L"$E_{0}(\alpha)$",
    xlab = L"\alpha", framestyle = :box, xtickfontsize=10,ytickfontsize=10, gridlinewidth = 0.0,
    guidefontsize = 14, linewidth=1.5, markerstrokewidth=1.0, legendfontsize = 10, legend=:topright,
    foreground_color_legend = nothing, extra_kwargs=Dict(:subplot=>Dict(:legend_hfactor=>1.25)))
    plot!(collect(1.0:16.0), VMC_energy.(collect(1.0:16.0)), label = "Feynman")
    plot!(collect(1.0:16.0), -collect(1.0:16.0).-1.26*(collect(1.0:16.0)./10).^2, label = "Perturbation")
    #vspan!([11,16], linecolor = :grey, fillcolor = :grey, label = "With Swap", alpha = 0.2)


    #savefig("C:/Users/60125/OneDrive - Imperial College London/Y4 Masters Project/Graphics/energy_1-16")
end

begin
    energy_list_total = vcat(energy_list, energy_list2, energy_list4[2:5])
    variance_list_total = vcat(variance_list, variance_list2, variance_list4[2:5])
    plot_font = "Computer Modern"
    default(fontfamily=plot_font, dpi = 500)

    plot(collect(1.0:20.0),  energy_list_total, yerr=sqrt.(variance_list_total), label = "DiagMC", ylab = L"$E_{0}(\alpha)$",
    xlab = L"\alpha", framestyle = :box, xtickfontsize=10,ytickfontsize=10, gridlinewidth = 0.0,
    guidefontsize = 14, linewidth=1.5, markerstrokewidth=1.0, legendfontsize = 10, legend=:topright,
    foreground_color_legend = nothing, extra_kwargs=Dict(:subplot=>Dict(:legend_hfactor=>1.25)))
    plot!(collect(1.0:20.0), VMC_energy.(collect(1.0:20.0)), label = "Feynman")
    plot!(collect(1.0:20.0), -collect(1.0:20.0).-1.26*(collect(1.0:20.0)./10).^2, label = "Perturbation")
    #vspan!([11,16], linecolor = :grey, fillcolor = :grey, label = "With Swap", alpha = 0.2)


    #savefig("C:/Users/60125/OneDrive - Imperial College London/Y4 Masters Project/Graphics/energy_1-20")
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
        diagram, estimators, variance = simulate!(diagram, estimators, false)
        push!(diagram_list2, diagram)
        push!(estimators_list2, estimators)
        push!(variance_list2, variance)
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
    pertub(k) = k^2/2 - (sqrt(2)/k)*asin(k/sqrt(2))-0.0026
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
    plot!(collect(range(0, sqrt(2), length=100)), pertub.(collect(range(0, sqrt(2), length=100))), label = "Perturbation", linecolor=:green)
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


begin
    α=3.0
    max_τ=100.0; max_order=5000;
    n_loop=500; n_hist=100000; sample_freq=200;
    k_list = collect(range(0, 2.0, length=25))
    diagram_k3list = []
    estimators_k3list = []
    variance_k3list = []
    for k in k_list
        diagram = initialise_diagram(α, k, max_τ, max_order)
        estimators = Estimators_Record(max_τ, max_order, sample_freq, n_loop, n_hist)
        diagram, estimators, variance = simulate!(diagram, estimators, false)
        push!(diagram_k3list, diagram)
        push!(estimators_k3list, estimators)
        push!(variance_k3list, variance)
    end
end

begin
    energy_k3list = []
    mass_k3list = []
    variance_k3list = []
    for estimator in estimators_k3list
        push!(energy_k3list, estimator.energy_mean[end])
        push!(mass_k3list, estimator.mass_mean[end])
        push!(variance_k3list, std(estimator.energy_mean[100000:end])^2)
    end
end

begin
    α=3.0
    max_τ=100.0; max_order=5000;
    n_loop=1000; n_hist=100000; sample_freq=200;
    k_list = collect(range(2.0, 3.0, length=20))
    diagram_k3list2 = []
    estimators_k3list2 = []
    variance_k3list2 = []
    for k in k_list
        diagram = initialise_diagram(α, k, max_τ, max_order)
        estimators = Estimators_Record(max_τ, max_order, sample_freq, n_loop, n_hist)
        diagram, estimators, variance = simulate!(diagram, estimators, false)
        push!(diagram_k3list2, diagram)
        push!(estimators_k3list2, estimators)
        push!(variance_k3list2, variance)
    end
end

begin
    energy_k3list2 = []
    mass_k3list2 = []
    variance_k3list2 = []
    for estimator in estimators_k3list2
        push!(energy_k3list2, estimator.energy_mean[end])
        push!(mass_k3list2, estimator.mass_mean[end])
        push!(variance_k3list2, std(estimator.energy_mean[100000:end])^2)
    end
end

begin
    energy_klist_total = vcat(energy_k3list, energy_k3list2[2:18])
    k_list_total = vcat(collect(range(0, 2.0, length=25)), collect(range(2.0, 3.0, length=20))[2:18])
    variance_klist_total = vcat(variance_k3list, variance_k3list2[2:18])
    quadratic(t, p) = p[1].+p[2].*t.+p[3].*t.*t
    p0=[0.0,-1.0,1.0]
    w=1 ./sqrt.(variance_klist_total[1:8])
    fit = curve_fit(quadratic, k_list_total[1:8], energy_klist_total[1:8], w,p0)
    plot_font = "Computer Modern"
    default(fontfamily=plot_font, dpi = 500)

    plot(k_list_total,  energy_klist_total, yerr=sqrt.(variance_klist_total), label = L"$\alpha=3$", ylab = L"$E_{0}(k)$",
    xlab = L"k", framestyle = :box, xtickfontsize=10,ytickfontsize=10, gridlinewidth = 0.0,
    guidefontsize = 14, linewidth=2.0, markerstrokewidth=1.0, legendfontsize = 10, legend=:topleft,
    foreground_color_legend = nothing, extra_kwargs=Dict(:subplot=>Dict(:legend_hfactor=>1.25)))
    plot!(k_list_total, quadratic(k_list_total, fit.param), label = "Parabolic fit", linestyle=:dash, linecolor=:red)
    hline!([energy_klist_total[1]+1], linestyle=:dot, label = "Continuum edge", linecolor=:grey)
    #plot!(alpha_plot_list, VMC_plot_list, label = "Feynman")


    savefig("C:/Users/60125/OneDrive - Imperial College London/Y4 Masters Project/Graphics/dispersion_alpha_3")
end

begin
    energy_klist_total = vcat(energy_klist[1:18], energy_klist2[1:8])
    k_list_total = vcat(collect(range(0, 2.0, length=25))[1:18], collect(range(1.5, 1.9, length=10))[1:8])
    variance_klist_total = vcat(variance_klist[1:18], variance_klist2[1:8])
    energy_klist_total2 = vcat(energy_k2list, energy_k2list2[2:8])
    k_list_total2 = vcat(collect(range(0, 2.0, length=25)), collect(range(2.0, 2.5, length=10))[2:8])
    variance_klist_total2 = vcat(variance_k2list, variance_k2list2[2:8])
    energy_klist_total3 = vcat(energy_k3list, energy_k3list2[2:18])
    k_list_total3 = vcat(collect(range(0, 2.0, length=25)), collect(range(2.0, 3.0, length=20))[2:18])
    variance_klist_total3 = vcat(variance_k3list, variance_k3list2[2:18])
    plot_font = "Computer Modern"
    default(fontfamily=plot_font, dpi = 500)


    plot(k_list_total,  energy_klist_total.-energy_klist_total[1], yerr=sqrt.(variance_klist_total), label = L"$\alpha=1$", ylab = L"$E_{0}(k)$",
    xlab = L"k", framestyle = :box, xtickfontsize=10,ytickfontsize=10, gridlinewidth = 0.0,
    guidefontsize = 14, linewidth=1.5, markerstrokewidth=1.5, legendfontsize = 10, legend=:topleft,
    foreground_color_legend = nothing, extra_kwargs=Dict(:subplot=>Dict(:legend_hfactor=>1.25)), markercolor=:grey)
    
    plot!(k_list_total2,  energy_klist_total2.-energy_klist_total2[1], yerr=sqrt.(variance_klist_total2), label = L"$\alpha=2$", ylab = L"$E_{0}(k)$",
    xlab = L"k", framestyle = :box, xtickfontsize=10,ytickfontsize=10, gridlinewidth = 0.0,
    guidefontsize = 14, linewidth=1.5, markerstrokewidth=1.5, legendfontsize = 10, legend=:topleft,
    foreground_color_legend = nothing, extra_kwargs=Dict(:subplot=>Dict(:legend_hfactor=>1.25)),markercolor=:grey)

    plot!(k_list_total3,  energy_klist_total3.-energy_klist_total3[1], yerr=sqrt.(variance_klist_total3), label = L"$\alpha=3$", ylab = L"$E_{0}(k,\alpha)-E_{0}(0,\alpha)$",
    xlab = L"k", framestyle = :box, xtickfontsize=10,ytickfontsize=10, gridlinewidth = 0.0,
    guidefontsize = 14, linewidth=1.5, markerstrokewidth=1.5, legendfontsize = 10, legend=:bottomright,
    foreground_color_legend = nothing, extra_kwargs=Dict(:subplot=>Dict(:legend_hfactor=>1.25)),markercolor=:grey)
    hline!([1], linestyle=:dash, label = "Continuum edge", linecolor=:grey, linewidth=1.0)
    #plot!(alpha_plot_list, VMC_plot_list, label = "Feynman")


    #savefig("C:/Users/60125/OneDrive - Imperial College London/Y4 Masters Project/Graphics/dispersion_alpha_123")
end

begin
    α_list = collect(range(0.2, 5.0, length=25))
    loop_list = [1000]
    p=0.0
    max_τ=100.0; max_order=5000;
    n_hist=100000; sample_freq=200;
    diagram_listmass = []
    estimators_listmass = []
    variance_listmass = []
    for α in α_list
        n_loop = loop_list[1]
        diagram = initialise_diagram(α, p, max_τ, max_order)
        estimators = Estimators_Record(max_τ, max_order, sample_freq, n_loop, n_hist)
        diagram, estimators, variance = simulate!(diagram, estimators, false)
        push!(diagram_listmass, diagram)
        push!(estimators_listmass, estimators)
        push!(variance_listmass, variance)
    end
end

begin
    energy_listmass = []
    mass_listmass = []
    variance_listmass = []
    for estimator in estimators_listmass
        push!(energy_listmass, estimator.energy_mean[end])
        push!(mass_listmass, estimator.mass_mean[end])
        push!(variance_listmass, std(estimator.mass_mean[100000:end])^2)
    end
end


begin


    mass_list_total = mass_listmass
    variance_list_total = variance_listmass
    plot_font = "Computer Modern"
    default(fontfamily=plot_font, dpi = 500)

    plot(collect(range(0.2, 5.0, length=25)),  mass_list_total, yerr=sqrt.(variance_list_total), label = "DiagMC", ylab = L"$m_{*}(\alpha)$",
    xlab = L"\alpha", framestyle = :box, xtickfontsize=10,ytickfontsize=10, gridlinewidth = 0.0,
    guidefontsize = 14, linewidth=1.5, markerstrokewidth=1.0, legendfontsize = 10, legend=:topleft,
    foreground_color_legend = nothing, extra_kwargs=Dict(:subplot=>Dict(:legend_hfactor=>1.25)))
    plot!(collect(range(0.2, 5.0, length=25)), polaron_effective_mass.(collect(range(0.2, 5.0, length=25))), label = "Feynman")
    plot!(collect(range(0.2, 5.0, length=25)), 1 ./(1 .-collect(range(0.2, 5.0, length=25))./6), label = "Perturbation", linestyle=:dash)
    #vspan!([11,16], linecolor = :grey, fillcolor = :grey, label = "With Swap", alpha = 0.2)


    #savefig("C:/Users/60125/OneDrive - Imperial College London/Y4 Masters Project/Graphics/eff_mass")
end