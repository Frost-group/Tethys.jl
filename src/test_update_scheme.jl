include("Diagram.jl")
include("update.jl")
include("measure.jl")
using Random
using LsqFit
using JLD2
using Logging


begin
    num_mea=1; regime=Diff_more(); regime_2=Diff_2()
    p=0; max_τ=30; max_order=500; mass=1; μ=-2.2; ω=1; α=2
end

begin
    hist=Hist_Record(300,max_τ,500)
    diagram=Diagram(0, max_τ, max_order, mass, μ, ω, α)
end

begin
    green_record,green_func,variance=hist_measure_4!(diagram,hist,10000)#
    @info "End of sampling"
end

begin
    max_time=100
    plot(hist.time_points[1:max_time],sum(green_func[i,1:max_time] for i in 1:501),yscale=:log)
end

begin
    statis=hist.normalized_data
    bin_width=hist.time_points[2]-hist.time_points[1]
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
    statis=hist.normalized_data
    min_time=Int(div(5,bin_width,RoundUp))
    max_time=Int(div(15,bin_width,RoundUp))
    time_1=collect(min_time:max_time)*bin_width.-(bin_width/2)
    data=[sum(statis[:,i]) for i in min_time:max_time]

    model(t, p) = p[1] * exp.(-p[2] * t)
    p0=[sum(statis[:,1]),(-α-μ)]
    fit = curve_fit(model, time_1, data, p0)

    linear(t, p) = p[1].-p[2].*t
    p0=[0,(-α-1.26*(α/10)^2-μ)]
    fit = curve_fit(linear, time_1, log.(data), p0)

    plot(time_1,data,markershape=:circle,yaxis=:log)

    display(plot!(time_1,model(time_1, [exp(fit.param[1]),fit.param[2]]),markershape=:circle,yaxis=:log))
    println(fit.param)
    # println("error ",standard_errors(fit))
    ener=-α-1.26*(α/10)^2
    println((ener-μ))
    println()
    println(fit.param[2]+μ)
    println(ener)
end

begin
    bin_width=hist.time_points[2]-hist.time_points[1]
    min_time=Int(div(5,bin_width,RoundUp))
    max_time=Int(div(15,bin_width,RoundUp))
    time_points=hist.time_points[min_time:max_time]
    statis=sum(green_func[i,:] for i in 1:max_order+1)
    y=log.(statis)[min_time:max_time]
    w=1 ./(variance[min_time:max_time]./(statis[min_time:max_time]).^2)
    linear(t, p) = p[1].-p[2].*t
    p0=[0,(-α-1.26*(α/10)^2-μ)]
    fit = curve_fit(linear, time_points, y, w, p0)
    println(fit.param[2]+μ)
    println(standard_errors(fit)[2])
end

begin
    plot(time_points,y)
end

begin
    n_loop=20000
    num_samples=5
    α_list=collect(1:num_samples)*0.3
    μ_list=-α_list.-1.26*(α_list./10).^2 .-0.5

    num_mea=1; regime=Diff_more();
    p=0; max_τ=30; max_order=500; mass=1; ω=1;

    energy_record=[]
    E_error_record=[]
    z0_record=[]
    Z_error_record=[]


    linear(t, p) = p[1].-p[2].*t
    bin_width=max_τ/300
    min_time=Int(div(6,bin_width,RoundUp))
    max_time=Int(div(10,bin_width,RoundUp))

    for i in 1:num_samples
        α=α_list[i]
        μ=μ_list[i]
        hist=Hist_Record(300,max_τ,max_order)
        diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)
        green_record,green_func,variance=hist_measure!(diagram,hist,n_loop)#
        println("end:",i)

        time_points=hist.time_points[min_time:max_time]
        statis=green_record[n_loop]
        y=log.(statis)[min_time:max_time]
        w=1 ./variance[min_time:max_time]
        p0=[0,(-α-1.26*(α/10)^2-μ)]
        fit = curve_fit(linear, time_points, y, w, p0)
        errors=standard_errors(fit)
        append!(energy_record,fit.param[2]+μ)
        append!(E_error_record,errors[2])

        z0=exp(fit.param[1])
        append!(z0_record,z0)
        append!(Z_error_record,z0*errors[1])
    end
end

begin
    scatter(α_list,energy_record,yerr=E_error_record,
    title="Ground state energy against α",label="DiagMC",titlefontsize=10)
    plot!(α_list,-α_list.-1.26*(α_list./10).^2,xlabel="α",ylabel="Energy",label="VMC")
    # savefig("energy.png")
end

begin
    scatter(α_list,z0_record,yerr=Z_error_record,
    xlabel="α",ylabel="Z₀",title="Z₀ against α",titlefontsize=10)
    # savefig("Z_0.png")
end

#μ test
begin
    n_loop=12000
    num_samples=5
    α_list=collect(1:num_samples)*0.3
    μ_list=-α_list.-1.26*(α_list./10).^2 .-0.5
    α=1.5
    μ_list=-(1.7 .+collect(0:num_samples-1).*0.1)
    num_mea=1; regime=Diff_more();
    p=0; max_τ=30; max_order=500; mass=1; ω=1;

    energy_record=[]
    E_error_record=[]
    z0_record=[]
    Z_error_record=[]


    linear(t, p) = p[1].-p[2].*t
    bin_width=max_τ/300
    min_time=Int(div(5,bin_width,RoundUp))
    max_time=Int(div(15,bin_width,RoundUp))

    for i in 1:num_samples
        # α=α_list[i]
        μ=μ_list[i]
        hist=Hist_Record(300,max_τ,max_order)
        diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)
        green_record,green_func,variance=hist_measure_3!(diagram,hist,n_loop)#
        println("end:",i)

        time_points=hist.time_points[min_time:max_time]
        # statis=green_record[n_loop]
        # y=log.(statis)[min_time:max_time]
        # w=1 ./variance[min_time:max_time]

        statis=sum(green_func[i,:] for i in 1:max_order+1)
        y=log.(statis)[min_time:max_time]
        w=1 ./(variance[min_time:max_time]./(statis[min_time:max_time]).^2)

        p0=[0,(-α-1.26*(α/10)^2-μ)]
        fit = curve_fit(linear, time_points, y, w, p0)
        errors=standard_errors(fit)
        append!(energy_record,fit.param[2]+μ)
        append!(E_error_record,errors[2])

        z0=exp(fit.param[1])
        append!(z0_record,z0)
        append!(Z_error_record,z0*errors[1])
        plot(time_points,y)
        display(plot!(time_points,linear(time_points,fit.param),
        xlabel="τ",ylabel="log(green)",title="μ="*string(μ)))
        # n_loop+=500
    end
end

begin 
    plot(μ_list,energy_record,yerr=E_error_record,xlabel="μ",ylabel="Energy",title="α="*string(α),label="DiagMC")
    plot!(μ_list,ones(num_samples)*(-α-1.26*(α/10)^2),xlabel="μ",ylabel="Energy",label="VMC")
end


begin
    disc=Discrete_Record(30,5,15,max_order)
    diagram=Diagram(0, max_τ, max_order, mass, μ, ω, α)
end

begin
    green_record,green_func=discrete_measure!(diagram,disc)
    println("end")
end

begin
    max_time=140
    plot(disc.time_points,sum(green_func[i,1:max_time] for i in 1:501),yscale=:log)
end