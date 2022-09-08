include("Diagram.jl")
include("update.jl")
include("measure.jl")
using Random
using LsqFit
using JLD2


begin
    num_mea=1; regime=Diff_more(); regime_2=Diff_2()
    p=0; max_τ=30; max_order=500; mass=1; μ=-2.2; ω=1; α=2
end

begin
    bin_width=max_τ/300
    hist=Hist_Record(300,max_τ,500)
    diagram=Diagram(0, max_τ, max_order, mass, μ, ω, α)
end

begin
    green_record,green_func,variance=hist_measure!(diagram,hist)#
    println("end")
end

begin
    plot(hist.time_points[1:150],sum(green_func[i,1:150] for i in 1:501),yscale=:log)
end

begin
    statis=hist.normalized_data
    max_orders=3
    max_time=50
    time_1=collect(1:max_time)*bin_width.-(bin_width/2)
    plot(time_1,statis[1,1:max_time],label="order:"*string(0))#,yaxis=:log,ylim=(0,5000000))
    for order in 2:max_orders+1
        if order == max_orders+1
            display(plot!(time_1,statis[order,1:max_time],label="order:"*string(order-1)))
        else
            plot!(time_1,statis[order,1:max_time],label="order:"*string(order-1))
        end
    end
end

begin
    min_time=Int(div(5,bin_width,RoundUp))
    max_time=Int(div(10,bin_width,RoundUp))
    time_1=collect(min_time:max_time)*bin_width.-(bin_width/2)
    data=[sum(statis[:,i]) for i in min_time:max_time]
    model(t, p) = p[1] * exp.(-p[2] * t)
    p0=[sum(statis[:,1]),(-α-μ)]
    fit = curve_fit(model, time_1, data, p0)
    plot(time_1,data,markershape=:circle,yaxis=:log)
    display(plot!(time_1,model(time_1, fit.param),markershape=:circle,yaxis=:log))
    println(fit.param)
    # println("error ",standard_errors(fit))
    ener=-α-1.26*(α/10)^2
    println((ener-μ))
    println()
    println(fit.param[2]+μ)
    println(ener)
end