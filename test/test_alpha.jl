using .Tethys
using Random
using LsqFit
using Plots
using JLD2

begin
    scanned_α=[]
    energy_record=[]
    E_error_record=[]
    z0_record=[]
    Z_error_record=[]
end

begin
    n_loop=10
    num_samples=30
    α_list=collect(26:26)*0.25.+0.1#num_samples-1 21 22
    μ_list=-α_list.-1.26*(α_list./10).^2 .-0.8
    # α=1.5
    # μ_list=-(1.7 .+collect(0:num_samples-1).*0.1)
    num_mea=1; regime=Diff_more();
    p=0; max_τ=30; max_order=500; mass=1; ω=1;


    linear(t, p) = p[1].-p[2].*t
    bin_width=max_τ/300
    min_time=Int(div(5,bin_width,RoundUp))
    max_time=Int(div(12,bin_width,RoundUp))

    for i in 1:num_samples
        α=α_list[i]
        μ=μ_list[i]
        # if i==22
        #     μ-=0.1
        # end
        # if i>22
        #     μ-=0.2
        # end
        hist=Hist_Record(300,max_τ,max_order)
        diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)
        diagram,hist,green_record,zero_record,green_func,variance=hist_measure!(diagram,hist,"E://data",true,n_loop)#
        println("end:",i)

        time_points=hist.time_points[min_time:max_time]

        statis=sum(green_func[i,:] for i in 1:max_order+1)
        y=log.(statis)[min_time:max_time]
        w=1 ./(variance[min_time:max_time]./(statis[min_time:max_time]).^2)

        p0=[0,(-α-1.26*(α/10)^2-μ)]
        fit = curve_fit(linear, time_points, y, w, p0)
        errors=standard_errors(fit)

        append!(scanned_α,α)
        append!(energy_record,fit.param[2]+μ)
        append!(E_error_record,errors[2])

        z0=exp(fit.param[1])
        append!(z0_record,z0)
        append!(Z_error_record,z0*errors[1])
        plot(time_points,y)
        display(plot!(time_points,linear(time_points,fit.param),
        xlabel="τ",ylabel="log(green)",title="α="*string(α)*",k="*string(0)*",μ="*string(μ)))
        # n_loop+=500
        println("energy:",fit.param[2]+μ)
        println("perturb:",-α-1.26*(α/10)^2)
    end
end

begin 
    plot(scanned_α,energy_record,yerr=E_error_record,xlabel="α",ylabel="Energy",label="DiagMC")
    plot!(scanned_α,-scanned_α.-1.26*(scanned_α./10).^2,xlabel="α",ylabel="Energy",label="Pertub")
end

begin 
    plot(scanned_α,z0_record,yerr=Z_error_record,xlabel="α",ylabel="Z_0",label="DiagMC")
end

# begin
#     n_loop=150000
#     num_samples=1
#     α_list=collect(1:1+num_samples)*0.8
#     μ_list=-α_list.-1.26*(α_list./10).^2 .-0.3

#     num_mea=1; regime=Diff_more();
#     p=0; max_τ=30; max_order=500; mass=1; ω=1;


#     linear(t, p) = p[1].-p[2].*t
#     bin_width=max_τ/300
#     min_time=Int(div(5,bin_width,RoundUp))
#     max_time=Int(div(12,bin_width,RoundUp))

#     for i in 1:num_samples
#         α=α_list[i]
#         μ=μ_list[i]
#         hist=Hist_Record(300,max_τ,max_order)
#         diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)
#         diagram,hist,green_record,zero_record,green_func,variance=hist_measure!(diagram,hist,"E://data",n_loop)#
#         println(sum(hist.unnormalized_data))
#         println("end:",i)

#         time_points=hist.time_points[min_time:max_time]

#         statis=sum(green_func[i,:] for i in 1:max_order+1)
#         y=log.(statis)[min_time:max_time]
#         w=1 ./(variance[min_time:max_time]./(statis[min_time:max_time]).^2)

#         p0=[0,(-α-1.26*(α/10)^2-μ)]
#         fit = curve_fit(linear, time_points, y, w, p0)
#         errors=standard_errors(fit)
#         append!(energy_record,fit.param[2]+μ)
#         append!(E_error_record,errors[2])

#         z0=exp(fit.param[1])
#         append!(z0_record,z0)
#         append!(Z_error_record,z0*errors[1])
#         plot(time_points,y)
#         display(plot!(time_points,linear(time_points,fit.param),
#         xlabel="τ",ylabel="log(green)",title="μ="*string(μ)))
#         # n_loop+=500
#         println("energy:",fit.param[2]+μ)
#         println("perturb:",-α-1.26*(α/10)^2)
#     end
# end

# begin
#     n_loop=150000
#     num_samples=1
#     α_list=collect(7:6+num_samples)*0.8
#     μ_list=-α_list.-1.26*(α_list./10).^2 .-0.4

#     num_mea=1; regime=Diff_more();
#     p=0; max_τ=30; max_order=500; mass=1; ω=1;


#     linear(t, p) = p[1].-p[2].*t
#     bin_width=max_τ/300
#     min_time=Int(div(5,bin_width,RoundUp))
#     max_time=Int(div(10,bin_width,RoundUp))

#     for i in 1:num_samples
#         α=α_list[i]
#         μ=μ_list[i]
#         hist=Hist_Record(300,max_τ,max_order)
#         diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)
#         diagram,hist,green_record,zero_record,green_func,variance=hist_measure!(diagram,hist,"E://data",n_loop)#
#         println("end:",i)

#         time_points=hist.time_points[min_time:max_time]

#         statis=sum(green_func[i,:] for i in 1:max_order+1)
#         y=log.(statis)[min_time:max_time]
#         w=1 ./(variance[min_time:max_time]./(statis[min_time:max_time]).^2)

#         p0=[0,(-α-1.26*(α/10)^2-μ)]
#         fit = curve_fit(linear, time_points, y, w, p0)
#         errors=standard_errors(fit)
#         append!(energy_record,fit.param[2]+μ)
#         append!(E_error_record,errors[2])

#         z0=exp(fit.param[1])
#         append!(z0_record,z0)
#         append!(Z_error_record,z0*errors[1])
#         plot(time_points,y)
#         display(plot!(time_points,linear(time_points,fit.param),
#         xlabel="τ",ylabel="log(green)",title="μ="*string(μ)))
#         # n_loop+=500
#         println("energy:",fit.param[2]+μ)
#         println("perturb:",-α-1.26*(α/10)^2)
#     end
# end

# begin 
#     α_list=[0.8,1.6,2.4,3.2,4.0,4.0,4.8,5.6,4.0,5.6,5.6,5.6,5.6]
#     plot(α_list,energy_record,yerr=E_error_record,xlabel="α",ylabel="Energy",label="DiagMC")
#     plot!(α_list,-α_list.-1.26*(α_list./10).^2,xlabel="α",ylabel="Energy",label="Pertub")
# end

# begin
#     order=[1,2,3,4,9,7,13]#,5,6]
#     α_list_1=[]
#     energy_record_1=[]
#     E_error_record_1=[]
#     z0_record_1=[]
#     Z_error_record_1=[]
#     for i in order
#         append!(α_list_1,α_list[i])
#         append!(energy_record_1,energy_record[i])
#         append!(E_error_record_1,E_error_record[i])
#         append!(z0_record_1,z0_record[i])
#         append!(Z_error_record_1,Z_error_record[i])
#     end
# end

# begin 
#     plot(α_list_1,energy_record_1,yerr=E_error_record_1,xlabel="α",ylabel="Energy",label="DiagMC")
#     plot!(α_list_1,-α_list_1.-1.26*(α_list_1./10).^2,xlabel="α",ylabel="Energy",label="Pertub")
# end

# begin 
#     plot(α_list_1,z0_record_1,yerr=Z_error_record_1,xlabel="α",ylabel="Z_0",label="DiagMC")
# end