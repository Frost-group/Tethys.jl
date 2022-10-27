include("Diagram.jl")
include("update.jl")
include("measure.jl")
using Random
using LsqFit
using JLD2
using Logging
using Base.Threads

begin
    num_mea=1; regime=Diff_more(); regime_2=Diff_2()
    p=0; max_τ=30; max_order=500; mass=1; μ=-2.2; ω=1; α=2
end

begin
    hist=Hist_Record(300,max_τ,500)
    diagram=Diagram(0, max_τ, max_order, mass, μ, ω, α)
end

begin
    n_loop=10000
    p_samples = 5
    p_list=collect(0:p_samples)*(1/p_samples)
    α=1
    μ=-1.2
    num_mea=1; regime=Diff_more();
    p=0; max_τ=30; max_order=500; mass=1; ω=1;

    energy_record = []
    E_error_record = []
    z0_record = []
    Z_error_record = []
    order_correction = []


    linear(t, p) = p[1].-p[2].*t
    bin_width=max_τ/300
    min_time=Int(div(5,bin_width,RoundUp))
    max_time=Int(div(12,bin_width,RoundUp))

    @threads for i in 1:p_samples+1
        p=p_list[i]
        hist=Hist_Record(300,max_τ,max_order)
        diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)
        green_record,green_func,variance=hist_measure_4!(diagram,hist,n_loop)#
        @info "End of loop" i

        time_points=hist.time_points[min_time:max_time]

        statis=sum(green_func[i,:] for i in 1:max_order+1)
        y=log.(statis)[min_time:max_time]
        w=1 ./(variance[min_time:max_time]./(statis[min_time:max_time]).^2)

        p0=[0,(-α-1.26*(α/10)^2-μ)]
        fit = curve_fit(linear, time_points, y, w, p0)
        errors=standard_errors(fit)
        append!(energy_record, fit.param[2]+μ)
        append!(E_error_record, errors[2])
        append!(order_correction, i)

        z0=exp(fit.param[1])
        append!(z0_record,z0)
        append!(Z_error_record,z0*errors[1])
        plot(time_points,y)
        display(plot!(time_points,linear(time_points,fit.param),
        xlabel="τ",ylabel="log(green)",title="p="*string(p)))
    end
end

begin
    quadratic(t, p) = p[1].*t.*t .+ p[2].*t .+ p[3]
    p0=[1,0.1,-1]
    fit = curve_fit(quadratic, p_list, energy_record[sortperm(order_correction)], p0)
    plot(p_list,energy_record[sortperm(order_correction)],yerr=E_error_record,xlabel="p",ylabel="Energy",title="Dispersion",label="DiagMC")
    display(plot!(p_list,linear(p_list,fit.param)))
end

