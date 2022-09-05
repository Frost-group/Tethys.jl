include("Diagram.jl")
include("update.jl")
using Random
using LsqFit
using JLD2

begin
    num_mea=1; regime=Diff_more(); regime_2=Diff_2()
    p_max=10; max_τ=30; max_order=500; mass=1; μ=-2.2; ω=1; α=2
end

begin
    statis=zeros(max_order+1, 300)
    bin_width=max_τ/300
    diagram_a=Diagram(0, max_τ, max_order, mass, μ, ω, α)
end

begin
    println("begin")
    # Random.seed!(1234)
    # record=[]
    loop=10000
    for j in 1:5000
        # loop_record=[]
        println(j)
        for i in 1:loop
            q=rand()
            # extend!(diagram_a)
            if diagram_a.order == 0
                if q<0.5
                    diagram_a.p_rem=0.2
                    insert_arc!(diagram_a,regime)
                    diagram_a.p_rem=0.5
                else
                    extend!(diagram_a)
                end
            else
                if q<0.2
                    insert_arc!(diagram_a,regime)
                elseif q<0.4
                    if diagram_a.order != 1
                        remove_arc!(diagram_a,regime) 
                    end       
                elseif q<0.6
                    swap_arc!(diagram_a)
                elseif q<0.8
                    extend!(diagram_a)
                else
                    if diagram_a.order == 1
                        diagram_a.p_rem=0.2
                        remove_arc!(diagram_a,regime) 
                        diagram_a.p_rem=0.5
                    end       
                end
            end
            statis[diagram_a.order+1,Int(div(diagram_a.τ,bin_width,RoundUp))]+=1
        end
    end
end

#zeroth order green function against time
begin
    max_time=60
    time_1=collect(1:max_time)*bin_width.-(bin_width/2)
    plot(time_1,statis[1,1:max_time],yaxis=:log,label="measured")
    display(plot!(time_1,statis[1,1]*exp.(μ.*time_1),label="theoretical", 
            title = "Zero order Green histogram μ="*string(μ),titlefontsize=12))
    println(-(log(statis[1,1])- log(statis[1,20]))/(bin_width*19))
    println(-(log(statis[1,75])- log(statis[1,150]))/(bin_width*74))
end

#green funcitons of different orders against time
begin
    max_orders=3
    max_time=50
    time_1=collect(1:max_time)*bin_width.-(bin_width/2)
    plot(time_1,statis[1,1:max_time],label="order:"*string(0))#,yaxis=:log,ylim=(0,5000000))
    # plot!(time_1,ones(max_time)*600)
    for order in 2:max_orders+1
        if order == max_orders+1
            display(plot!(time_1,statis[order,1:max_time],label="order:"*string(order-1)))
        else
            plot!(time_1,statis[order,1:max_time],label="order:"*string(order-1))
        end
    end
end

#number of data in different orders
begin
    max_check=6
    data=[sum(statis[i,:]) for i in 1:max_check+1]
    plot(collect(0:max_check),data, markershape=:circle,yaxis=:log,xticks=0:max_check)
end

#total green function agaisnt time
begin
    max_time=129#218
    time_1=collect(1:max_time)*bin_width.-(bin_width/2)
    data=[sum(statis[:,i]) for i in 1:max_time]
    display(plot(time_1,data,markershape=:circle,yaxis=:log))
    println(-(log(data[80])- log(data[60]))/(bin_width*59)+μ)
end

#fitting for ground state energy
begin
    min_time=Int(div(5,bin_width,RoundUp))
    max_time=Int(div(8,bin_width,RoundUp))
    time_1=collect(min_time:max_time)*bin_width.-(bin_width/2)
    data=[sum(statis[:,i]) for i in min_time:max_time]
    model(t, p) = p[1] * exp.(-p[2] * t)
    p0=[sum(statis[:,1]),(α-μ)]
    fit = curve_fit(model, time_1, data, p0)
    println(fit.param)
    println("error ",standard_errors(fit))
    ener=-α-1.26*(α/10)^2
    println((ener-μ))
    println()
    println(fit.param[2]+μ)
    println(ener)
end