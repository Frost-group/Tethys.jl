include("Diagram.jl")
include("update.jl")
include("measure.jl")

using Logging
using BenchmarkTools
using LsqFit

begin
    n_loop=100
    num_mea=1; regime=Diff_more();
    p=0; max_τ=30; max_order=500; mass=1; ω=1;


    linear(t, p) = p[1].-p[2].*t
    bin_width=max_τ/300
    min_time=Int(div(5,bin_width,RoundUp))
    max_time=Int(div(12,bin_width,RoundUp))

    α=20
    μ=-47
    n_hist = 100000
    directory = "D://data"
    store_data = false
    hist=Hist_Record(300,max_τ,max_order)
    diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)
    @info string(diagram.τ)
    diagram,hist,green_record,zero_record,green_func,variance=hist_measure!(diagram,hist,directory,n_loop,store_data,n_hist)

    #println(diagram.order)
    #println(diagram.τ)
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

begin
    plot(time_points,y)
end

begin
    plot(hist.time_points,log.(statis))
end

begin
    plot([0:1:1000;],log.(sum(green_func,dims=2)))
end

begin
    plot([0:1:500;],log.(sum(normalized_data,dims=2)))
end

begin
    linear(t, p) = p[1].-p[2].*t
    bin_width=max_τ/300
    min_time=Int(div(3,bin_width,RoundUp))
    max_time=Int(div(15,bin_width,RoundUp))
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

begin
    n_loop=1000
    num_mea=1; regime=Diff_more();
    p_list = [0.0:0.2:3.0;]
    max_τ=30; max_order=500; mass=1; ω=1;


    linear(t, p) = p[1].-p[2].*t
    bin_width=max_τ/300
    min_time=Int(div(5,bin_width,RoundUp))
    max_time=Int(div(12,bin_width,RoundUp))

    α=1
    μ=-1.2
    n_hist = 100000
    directory = "D://data"
    store_data = false
    for p in p_list
        timings = []
        μ_list = -α.-1.26*(α./10).^2 .+[-2.0:0.1:-0.1;]
        for μ in μ_list
            @info "Loading μ = "*string(μ)
            hist=Hist_Record(300,max_τ,max_order)
            diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)
            stats = @timed hist_measure!(diagram,hist,directory,n_loop,store_data,n_hist)
            append!(timings, stats.time/n_loop)

            @info "Time per sweep w/o data collection" stats.time/n_loop
            @info "Memory allocated" stats.bytes
        end
        display(plot(μ_list,timings,xlabel="μ",ylabel="Seconds per sweep",label="k = "*string(p)))
    end

end
