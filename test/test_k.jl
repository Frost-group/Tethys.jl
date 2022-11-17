using Tethys
using Random
using LsqFit
using Plots
using JLD2

begin
    scanned_k=[]
    energyk_record=[]
    Ek_error_record=[]
    zk_record=[]
    Zk_error_record=[]
end

begin
    n_loop=20000
    num_samples=20
    k_list=collect(15:num_samples-1)*0.15
    α=1
    μ_list=(k_list.^2)./(2*(1+α/6)) .+(-α-1.26*(α/10)^2-0.3)
    # α=1.5
    # μ_list=-(1.7 .+collect(0:num_samples-1).*0.1)
    num_mea=1; regime=Diff_more();
    p=0; max_τ=30; max_order=500; mass=1; ω=1;


    linear(t, p) = p[1].-p[2].*t
    bin_width=max_τ/300
    min_time=Int(div(5,bin_width,RoundUp))
    max_time=Int(div(10,bin_width,RoundUp))

    for i in 1:num_samples
        p=k_list[i]
        μ=-0.5
        hist=Hist_Record(300,max_τ,max_order)
        diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)
        diagram,hist,green_record,zero_record,green_func,variance=hist_measure!(diagram,hist,"E://data",true,n_loop)#
        println("end:",i)

        time_points=hist.time_points[min_time:max_time]

        statis=sum(green_func[i,:] for i in 1:max_order+1)
        y=log.(statis)[min_time:max_time]
        w=1 ./(variance[min_time:max_time]./(statis[min_time:max_time]).^2)

        p0=[0,((p^2)/(2*(1+α/6))-α-1.26*(α/10)^2-μ)]
        fit = curve_fit(linear, time_points, y, w, p0)
        errors=standard_errors(fit)

        append!(scanned_k,p)
        append!(energyk_record,fit.param[2]+μ)
        append!(Ek_error_record,errors[2])

        z0=exp(fit.param[1])
        append!(zk_record,z0)
        append!(Zk_error_record,z0*errors[1])
        plot(time_points,y)
        display(plot!(time_points,linear(time_points,fit.param),
        xlabel="τ",ylabel="log(green)",title="α="*string(α)*",k="*string(p)*",μ="*string(μ)))
        # n_loop+=500
        println("energy:",fit.param[2]+μ)
        println("perturb:",(p^2)/(2*(1+α/6))-α-1.26*(α/10)^2)

        # if i+1<=num_samples && fit.param[2]+μ<μ_list[i+1]
        #     μ_list[i+1]=fit.param[2]+μ
        # end

    end
end

begin 
    plot(scanned_k,energyk_record,yerr=Ek_error_record,xlabel="k",ylabel="Energy",label="DiagMC")
    plot!(scanned_k,(scanned_k.^2)./(2*(1+α/6)) .+(-α-1.26*(α/10)^2),xlabel="k",ylabel="Energy",label="Parabolic",
            title="α="*string(α)*" energy dispersion graph")
end

begin 
    plot(scanned_k,zk_record,yerr=Zk_error_record,xlabel="k",ylabel="Z_k",label="DiagMC",
        title="α="*string(α)*" quasi weight graph")
end

# begin
#     n_loop=150000
#     num_samples=1
#     k_list=collect(7:7)*0.5
#     k_list=[3.3]
#     α=1
#     μ_list=(k_list.^2)./(2*(1+α/6)) .+(-α-1.26*(α/10)^2-0.4)
#     # α=1.5
#     μ_list=[0.9]
#     # μ_list=-(1.7 .+collect(0:num_samples-1).*0.1)
#     num_mea=1; regime=Diff_more();
#     p=0; max_τ=30; max_order=500; mass=1; ω=1;

#     linear(t, p) = p[1].-p[2].*t
#     bin_width=max_τ/300
#     min_time=Int(div(5,bin_width,RoundUp))
#     max_time=Int(div(12,bin_width,RoundUp))

#     for i in 1:num_samples
#         p=k_list[i]
#         μ=μ_list[i]
#         hist=Hist_Record(300,max_τ,max_order)
#         diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)
#         green_record,green_func,variance=hist_measure_4!(diagram,hist,n_loop)#
#         println("end:",i)

#         time_points=hist.time_points[min_time:max_time]

#         statis=sum(green_func[i,:] for i in 1:max_order+1)
#         y=log.(statis)[min_time:max_time]
#         w=1 ./(variance[min_time:max_time]./(statis[min_time:max_time]).^2)
#         display(plot(time_points,y))

#         p0=[0,((p^2)/(2*(1+α/6))-α-1.26*(α/10)^2-μ)]
#         fit = curve_fit(linear, time_points, y, w, p0)
#         errors=standard_errors(fit)
#         append!(energyk_record,fit.param[2]+μ)
#         append!(Ek_error_record,errors[2])

#         z0=exp(fit.param[1])
#         append!(zk_record,z0)
#         append!(Zk_error_record,z0*errors[1])
#         plot(time_points,y)
#         display(plot!(time_points,linear(time_points,fit.param),
#         xlabel="τ",ylabel="log(green)",title="μ="*string(μ)))
#         # n_loop+=500
#         println("energy:",fit.param[2]+μ)
#         println("perturb:",(p^2)/(2*(1+α/6))-α-1.26*(α/10)^2)
        
#         if i+1<=num_samples && fit.param[2]+μ<μ_list[i+1]
#             μ_list[i+1]=fit.param[2]+μ
#         end

#     end
# end

# begin 
#     k_list=[0,0.5,1,1.5,2,2.5,3,2,2.5,3,
#             2.5,2.5,3,3,2,2.5,0.5,3.3,3.3,
#             3.3,3.3,3.3,3.3,3.3]
#     order=[1,17,3,4,15,16,14,22]#,18] #,5,6]
#     k_list_1=[]
#     energyk_record_1=[]
#     Ek_error_record_1=[]
#     zk_record_1=[]
#     Zk_error_record_1=[]
#     for i in order
#         append!(k_list_1,k_list[i])
#         append!(energyk_record_1,energyk_record[i])
#         append!(Ek_error_record_1,Ek_error_record[i])
#         append!(zk_record_1,zk_record[i])
#         append!(Zk_error_record_1,Zk_error_record[i])
#     end
# end

# begin 
#     plot(k_list_1,energyk_record_1,yerr=Ek_error_record_1,xlabel="k",ylabel="Energy",label="DiagMC")
#     plot!(k_list_1,(k_list_1.^2)./(2*(1+α/6)) .+(-α-1.26*(α/10)^2),xlabel="k",ylabel="Energy",label="Parabolic")
# end

# begin 
#     plot(k_list_1,zk_record_1,yerr=Zk_error_record_1,xlabel="k",ylabel="Z_k",label="DiagMC")
# end

begin
    using DataFrames
    a=[[1,2],[3,4]]
    # a=[hcat(i) for i in a]
    a=reshape(a, length(a), 1)
    a=DataFrame(a, :auto)
end