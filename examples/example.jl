include("../src/new_measure.jl")
using Random
using LsqFit
using JLD2
using FFTW
using Logging
using Statistics
using LaTeXStrings

"""
Code block to initialise and run a simulation.
Relevant parameters:
Electron-phonon coupling = 1
Diagram momentum = 0
Number of MC sweeps = 500
Number of diagram updates per sweep = 100000
Sampling frequency of measurements = 100

Swap scheme is set to false.
"""
begin
    α=1.0; p=0.0
    max_τ=1000.0; max_order=5000;
    n_loop=500; n_hist=100000; sample_freq=100;
    diagram = initialise_diagram(α, p, max_τ, max_order)
    estimators = Estimators_Record(max_τ, max_order, sample_freq, n_loop, n_hist)
    diagram, estimators, variance = simulate!(diagram, estimators,false)
end

"""
Plotting the mean energy evolution of a diagram with a coupling strength of 1.
"""
begin
    plot_font = "Computer Modern"
    default(fontfamily=plot_font, dpi = 500)

    plot(estimators.energy_mean[10:end], label = "DiagMC", ylab = L"$\langle E_{0}(\alpha=1)\rangle$",
    xlab = "Samples", framestyle = :box, xtickfontsize=8,ytickfontsize=9, color=:red,
    xguidefontsize = 12, linewidth=1.5, markerstrokewidth=1.0, legendfontsize = 11, legend=:topright,
    foreground_color_legend = nothing, extra_kwargs=Dict(:subplot=>Dict(:legend_hfactor=>1.25)),
    yguidefontsize=14)
    hline!([VMC_energy(1)], linestyle=:dash, label = "Feynman", linecolor=:black, linewidth=1.0)
end