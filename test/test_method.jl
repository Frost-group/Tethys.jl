#using .Tethys
include("../src/Diagram.jl")
include("../src/update.jl")
include("../src/zero_update.jl")
include("../src/measure.jl")
using Random
using LsqFit
using JLD2

begin
    n_loop=200000
    num_samples=20
    n_hist=100000
    α=1#15#20#20#5#20
    μ=0#-48#-47.1#-6.3#-47#3.5#9.2
    num_mea=1; regime=Diff_more();
    p=0; max_τ=30; max_order=2000; mass=1; ω=1;
    diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)
    set_τ!(diagram,20.0)
end

begin
    n_loop=400000
    p_ins=0.2;p_rem=0.2;p_from_0=1;
    real_normalized=[p_ins,p_rem]
    real_normalized/=sum(real_normalized)
    fake_normalized=[p_from_0]
    fake_normalized/=sum(fake_normalized)
    real_cumsum=cumsum(real_normalized)
    fake_cumsum=cumsum(fake_normalized)
    diagram.p_ins=real_normalized[1]
    diagram.p_rem=real_normalized[2]
    num_bins=300
    unnormalized_data=zeros(num_bins)
    bin_width=max_τ/num_bins
    weight_box = zeros(max_order)
    order_box= zeros(max_order)
    energy_record=[]
    # n_loop=1
    # n_hist=5000
end


begin
    Random.seed!(132432333)
    num_samples=1
    dia_order=diagram.order
    m=diagram.mass
    μ=diagram.μ
    ω=diagram.ω
    α=diagram.α
    α_squared=2pi*α*sqrt(2)
    result=true
    println("begin")
    for j in 1:n_loop
        println("loop.number:",j)
        for i in 1:n_hist
            q=rand()
            # if !result
            #     swap_arc!(diagram)
            # end
            if dia_order == 0
                diagram.p_ins=fake_normalized[1]
                result=insert_arc!(diagram,dia_order,m,μ,ω,α_squared)
                diagram.p_ins=real_normalized[1]
            elseif  dia_order == 1
                if q<real_cumsum[1]
                    result=insert_arc!(diagram,dia_order,m,μ,ω,α_squared)
                else
                    diagram.p_ins=fake_normalized[1]
                    result=remove_arc!(diagram,dia_order,m,μ,ω,α_squared)
                    diagram.p_ins=real_normalized[1]
                end
            else
                if q<real_cumsum[1]
                    result=insert_arc!(diagram,dia_order,m,μ,ω,α_squared)
                else
                    result=remove_arc!(diagram,dia_order,m,μ,ω,α_squared)
                end
            end
            # println()
            # check_k(diagram)
            # println()
            dia_order=diagram.order
            resample_arc!(diagram,dia_order,m,μ,ω,α_squared)
            E_value=energy(diagram)
            append!(energy_record,E_value)
            # swap_arc!(diagram)
            # println(diagram.arc_box)
            # check_timeorder(diagram)
            # println()
            # println("new")
            # println(diagram.component)
            # println(diagram.line_box)
            # println(diagram.arc_box)
            # println(diagram.end_arc_box)
            # println()
            # if resample_arc!(diagram,dia_order,m,μ,ω,α_squared)
            # # println(diagram.line_box)
            # # println(diagram.arc_box)
            #     println("yes")
            #     println("dis_t1")
            #     println(diagram.dispersion)
            #     println("dis_t2")
            #     println(total_dis_check(diagram))
            # else
            #     println("no")
            #     println("dis_t1")
            #     println(diagram.dispersion)
            #     println("dis_t2")
            #     println(total_dis_check(diagram))
            # end
            # println(diagram.line_box)
            # println(diagram.arc_box)
            # println(diagram.end_arc_box)
            # println()
            # check_timeorder(diagram)
            
            # if !result
            #     swap_arc!(diagram)
            # end
            
            # swap_arc!(diagram)
            
            # update_arcp!(diagram,dia_order,m,μ)
            # update_arcp!(diagram,dia_order,m,μ)
            # update_arcp!(diagram,dia_order,m,μ)


            # dis=-diagram.dispersion
            # println(dis)
            # if abs(2*dia_order-dis)<sqrt(dis)
            #     scale!(diagram,dia_order,m,μ,ω,num_samples)
            #     τ=diagram.record_τ[1]
            #     unnormalized_data[Int(div(τ,bin_width,RoundUp))]+=1
            # end

            # scale!(diagram,dia_order,m,μ,ω,num_samples)
            # # println(-diagram.dispersion/diagram.τ)
            # for i in 1:num_samples
            #     τ=diagram.record_τ[i]
            #     if τ<max_τ
            #         unnormalized_data[Int(div(τ,bin_width,RoundUp))]+=1#.0/pdf(Gamma(2*diagram.order+1,-7.5/diagram.dispersion),7.5)
            #     end
            # end
            component=diagram.component
            weight_box[component+1]+=1
            order_box[dia_order+1]+=1
            
            # E_value=energy(diagram)
            # append!(energy_record,E_value)
            # if mod(i,200) == 0
            #     E_value=energy(diagram)
            #     append!(energy_record,E_value)
            # end
        end
    end
end

begin
    plot(weight_box/sum(weight_box),xlims = (0,10),title="α="*string(α))
end

begin
    plot(order_box/sum(order_box),xlims = (0,100),title="α="*string(α))
    # plot!(range(0,100),pdf(Poisson(-diagram.dispersion/1.99),range(0,00)))
end

begin
    histogram(energy_record)#,xlims = (-5.0,-2.5))
    # println(mean(energy_record))
    # println(std(energy_record))
end

begin
    time=check_timeorder(diagram)
end

begin
    for i in 2:length(time)
        if time[i][1]==time[i-1][2]
            println(true)
        else
            println(false)
        end
    end
end

begin
    max_order=2000
    diagram.max_order=max_order
    n_loop=200000
    num_samples=20
    n_hist=100000
    unnormalized_data=zeros(num_bins)
    bin_width=max_τ/num_bins
    weight_box = zeros(max_order)
    # n_loop=200000000
    # n_hist=200
end

begin
    # Random.seed!(1353)
    num_samples=1
    dia_order=diagram.order
    m=diagram.mass
    μ=diagram.μ
    ω=diagram.ω
    α=diagram.α
    α_squared=2pi*α*sqrt(2)

    println("begin")
    for j in 1:n_loop
        println("loop.number:",j)
        for i in 1:n_hist
            q=rand()
            if dia_order == 0
                diagram.p_ins=fake_normalized[1]
                zero_insert_arc!(diagram,dia_order,m,μ,ω,α_squared)
                diagram.p_ins=real_normalized[1]
            elseif  dia_order == 1
                if q<real_cumsum[1]
                    zero_insert_arc!(diagram,dia_order,m,μ,ω,α_squared)
                else
                    diagram.p_ins=fake_normalized[1]
                    zero_remove_arc!(diagram,dia_order,m,μ,ω,α_squared)
                    diagram.p_ins=real_normalized[1]
                end
            else
                if q<real_cumsum[1]
                    zero_insert_arc!(diagram,dia_order,m,μ,ω,α_squared)
                else
                    zero_remove_arc!(diagram,dia_order,m,μ,ω,α_squared)
                end
            end
            swap_arc!(diagram)
            dia_order=diagram.order
            update_arcp!(diagram,dia_order,m,μ)

            # if mod(i,100) == 0
            #     dis=-diagram.dispersion
            #     if abs(2*dia_order-dis)<2*sqrt(dis)
            #         scale!(diagram,dia_order,m,μ,ω,num_samples)
            #         τ=diagram.record_τ[1]
            #         unnormalized_data[Int(div(τ,bin_width,RoundUp))]+=1
            #     end
            # end
            if mod(i,100) == 0
                scale!(diagram,dia_order,m,μ,ω,num_samples)
                τ=diagram.record_τ[1]
                unnormalized_data[Int(div(τ,bin_width,RoundUp))]+=1
            end

            component=diagram.component
            weight_box[component+1]+=1
            order_box[dia_order+1]+=1
        end
            # scale!(diagram,dia_order,m,μ,ω,num_samples)
            # # println(-diagram.dispersion/diagram.τ)
            # for i in 1:num_samples
            #     τ=diagram.record_τ[i]
            #     if τ<max_τ
            #         unnormalized_data[Int(div(τ,bin_width,RoundUp))]+=1#.0/pdf(Gamma(2*diagram.order+1,-7.5/diagram.dispersion),7.5)
            #     end
            # end
    end
end

begin
    time_points=collect(1:num_bins)*bin_width.-(bin_width/2)
    display(plot(time_points,log.(unnormalized_data),xlims = (0,10)))
    linear(t, p) = p[1].-p[2].*t
    min_time=Int(div(5.5,bin_width,RoundUp))
    max_time=Int(div(6,bin_width,RoundUp))
    p0=[0,(-α-1.26*(α/10)^2-μ)]
    y=log.(unnormalized_data)[min_time:max_time]
    time_points=time_points[min_time:max_time]
    fit = curve_fit(linear, time_points, y, p0)
end

begin
    min_time=Int(div(3,bin_width,RoundUp))
    max_time=Int(div(10,bin_width,RoundUp))
    time_points=collect(1:num_bins)*bin_width.-(bin_width/2)
    time_points=time_points[min_time:max_time]
    y=log.(unnormalized_data[min_time:max_time])
    second=first_dif(y,bin_width)
    plot(time_points[2:length(y)-1],second)
end

begin
    min_time=Int(div(4,bin_width,RoundUp))
    max_time=Int(div(5,bin_width,RoundUp))
    time_points=collect(1:num_bins)*bin_width.-(bin_width/2)
    time_points=time_points[min_time:max_time]
    y=unnormalized_data[min_time:max_time]
    display(plot(time_points,log.(y),xlims = (5,10)))
    linear(t, p) = p[1].-p[2].*t
    p0=[0,(-α-1.26*(α/10)^2-μ)]
    y=log.(y)
    fit = curve_fit(linear, time_points, y, p0)
end

begin
    plot(weight_box/sum(weight_box),xlims = (0,1000),title="α="*string(α))
end

begin
    set_μ!(diagram,-46.0)
end


function second_dif(data,bin_width)
    diff=[]
    for i in 2:length(data)-1
        append!(diff,(data[i-1]-2*data[i]+data[i+1])/bin_width^2)
    end
    return diff
end

function first_dif(data,bin_width)
    diff=[]
    for i in 2:length(data)-1
        append!(diff,(data[i+1]-data[i-1])/(2*bin_width))
    end
    return diff
end