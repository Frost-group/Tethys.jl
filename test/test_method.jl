#using .Tethys
include("../src/Diagram.jl")
include("../src/update.jl")
include("../src/zero_update.jl")
include("../src/measure.jl")
using Random
using LsqFit
using JLD2
using FFTW
using Logging
using Statistics
using LaTeXStrings

begin
    n_loop=200000
    num_samples=20
    n_hist=50000
    α=20.0#20#5#20
    μ=0#-48#-47.1#-6.3#-47#3.5#9.2
    num_mea=1; regime=Diff_more();
    p=0; max_τ=30; max_order=2000; mass=1; ω=1;
    diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)
    set_τ!(diagram,30.0)
end

begin
    p_record=[]
    record=[]
    order_se=[]
end

begin
    p=2
    diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)
    append!(p_record,p)
end


begin
    α=1#20#20.0
    diagram.α=α
    set_τ!(diagram,20.0)
    n_loop=200000
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
    c1=0
    r1=0
    c2=0
    r2=0
    c3=0
    r3=0
    # n_loop=1
    # n_hist=10#5000
end


begin
    # Random.seed!(132432333)
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
        order_box= zeros(max_order)
        for i in 1:n_hist
            q=rand()
            # if !result
            #     swap_arc!(diagram)
            # end
            if dia_order == 0
                c1+=1
                diagram.p_ins=fake_normalized[1]
                result=insert_arc!(diagram,dia_order,m,μ,ω,α_squared)
                diagram.p_ins=real_normalized[1]
                if result
                    r1+=1
                end
            elseif  dia_order == 1
                if q<real_cumsum[1]
                    c1+=1
                    result=insert_arc!(diagram,dia_order,m,μ,ω,α_squared)
                    if result
                        r1+=1
                    end
                else
                    c2+=1
                    diagram.p_ins=fake_normalized[1]
                    result=remove_arc!(diagram,dia_order,m,μ,ω,α_squared)
                    diagram.p_ins=real_normalized[1]
                    if result
                        r2+=1
                    end
                    #add = true
                end
            else
                if q<real_cumsum[1]
                    c1+=1
                    result=insert_arc!(diagram,dia_order,m,μ,ω,α_squared)
                    if result
                        r1+=1
                    end

                    # dia_order=diagram.order
                    # c1+=1
                    # result=insert_arc!(diagram,dia_order,m,μ,ω,α_squared)
                    # if result
                    #     r1+=1
                    # end
                else
                    c2+=1
                    result=remove_arc!(diagram,dia_order,m,μ,ω,α_squared)
                    if result
                        r2+=1
                    end

                    # dia_order=diagram.order
                    # c2+=1
                    # result=remove_arc!(diagram,dia_order,m,μ,ω,α_squared)
                    # if result
                    #     r2+=1
                    # end
                    #add = true
                end
            end
            # println()
            # check_k(diagram)
            # println()
            dia_order=diagram.order
            # update_arcp!(diagram,dia_order,m,μ)
            # c3+=1
            # result=resample_arc!(diagram,dia_order,m,μ,ω,α_squared)
            # if result
            #     r3+=1
            # end
            E_value=energy(diagram)
            append!(energy_record,E_value)
            swap_arc!(diagram)
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
            # println("dis_t1")
            # println(diagram.dispersion)
            # println("dis_t2")
            # println(total_dis_check(diagram))
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
        # push!(order_se,order_box)
    end
end

begin
    plot(weight_box/sum(weight_box),xlims = (0,5),title="α="*string(α))
    #push!(order_se,weight_box/sum(weight_box))
end

begin
    plot(order_box/sum(order_box),xlims = (0,200),title="α="*string(α))
    # plot!(range(0,100),pdf(Poisson(-diagram.dispersion/1.99),range(0,00)))
end

begin
    # order_box=order_se[1]
    # plot(order_box/sum(order_box),xlims = (0,1000),title="α="*string(α))
    p = plot()
    i=1
    for order_box in order_se
        p=plot!(range(0,length(order_box)-1),order_box,xlims = (0,8),title=L"α=1",xlabel=L"N",ylabel=L"Z_N",label=L"k=%$i"*string(p_record[i]),legend=:topright, dpi=300)
        i+=1
    end
    display(p)
    # savefig(p,"C://Users//wenzhaoren//Desktop//k_qua.png")
end

begin
    println(mean(energy_record[2000:end]))
    p=histogram(energy_record[2000:end],normalize=:probability,xlabel=L"E_0",ylabel="Frequency", legend=false,title=L"α=1",dpi=300)#,xlims = (-5.0,-2.5))
    # println(mean(energy_record))
    # println(std(energy_record))
    savefig(p,"C://Users//wenzhaoren//Desktop//freq.png")
end

begin
    push!(record,[α,mean(energy_record[2000:end])])
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
    for arc in diagram.arc_box
        println(norm(arc.q))
    end
    for arc in diagram.end_arc_box
        println(norm(arc.q))
    end
end

begin
    num = 0.0
    denom = 0.0
    ratio_list = []
    for val in arc_q_total_list
        if val == 0.0
            continue
        elseif val > rand()
            num += 1.0
            denom += 1.0
        else
            denom += 1.0
        end
        append!(ratio_list, num/denom)
    end
end

begin
    arc_q_mean = []
    for i in 1:length(arc_q_total_list)
        append!(arc_q_mean,mean(arc_q_total_list[1:i]))
    end
end

function auto_window(taus, c::Float64)
    m = collect(1:length(taus)) .< c * taus
    if sum(m) > 0
        return argmin(m)
    end
    return length(taus)
end

begin
    cutoff = 20000

    n = Int64(log2(nextpow(2,length(energy_record[cutoff:end]))))

    F = fft(energy_record[cutoff:end].-mean(energy_record[cutoff:end]))
    acf = real.(ifft(F.*conj.(F)))
    acf = acf/(4*n)
    acf = acf/acf[1]
    #freqs = fftshift(fftfreq(length(t), fs))

    #display(plot(collect(0: length(energy_record[cutoff:end])-1),acf))

    taus = 2.0 * cumsum(acf) .- 1.0

    window = auto_window(taus, 5.0)

end

function jackknife_energy(energy_samples, n_loop, n_hist, sample_freq)
    energy_array = deepcopy(energy_samples)

    samples_per_loop = Int64(n_loop*n_hist/sample_freq)
    
    energy_array = reshape(energy_array, (n_loop, div(length(energy_array), n_loop)))
    energy_array = [energy_array[:,i] for i in 1:size(energy_array,2)]
    energy_array = mean(energy_array, dims=1)[1]
    energy_array_sum = sum(energy_array)

    mean_value = mean(energy_array, dims=1)[1]
    
    variance_jk = 0.0
    for k in 1:n_loop
        jk_estimator = 0.0
        jk_estimator = energy_array_sum-energy_array[k]     
        jk_estimator *= 1/(n_loop-1)
        variance_jk += (jk_estimator - mean_value)^2

    end
    variance_jk *= (n_loop-1)/n_loop

    return variance_jk
end

function jackknife_energy_test(energy_samples, n_loop, n_hist, sample_freq)
    energy_array = deepcopy(energy_samples)
    energy_array_sum = sum(energy_array)
    mean_value = mean(energy_array, dims=1)[1]
    
    variance_jk = 0.0
    for k in 1:length(energy_array)
        jk_estimator = 0.0
        jk_estimator = energy_array_sum-energy_array[k]
        jk_estimator *= 1/(length(energy_array)-1)
        variance_jk += (jk_estimator - mean_value)^2

    end
    variance_jk *= (length(energy_array)-1)/length(energy_array)

    return variance_jk
end

function dispersion_relation(α, p, n_loop)
    n_hist=100000
    μ=0
    max_τ=30; max_order=2000; mass=1; ω=1;
    diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)
    p_ins=0.2;p_rem=0.2;p_from_0=1;
    real_normalized=[p_ins,p_rem]
    real_normalized/=sum(real_normalized)
    fake_normalized=[p_from_0]
    fake_normalized/=sum(fake_normalized)
    real_cumsum=cumsum(real_normalized)
    diagram.p_ins=real_normalized[1]
    diagram.p_rem=real_normalized[2]
    weight_box = zeros(max_order)
    order_box= zeros(max_order)
    sample_freq = 200
    energy_record=zeros(Int64(n_loop*n_hist/sample_freq))

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

            #if !result
            #    swap_arc!(diagram)
            #end
            #swap_arc!(diagram)
            dia_order=diagram.order
            component=diagram.component
            weight_box[component+1]+=1
            order_box[dia_order+1]+=1
            if mod(i,sample_freq) == 0
                E_value=energy(diagram)
                energy_record[Int64(((j-1)*n_hist+i)/sample_freq)] = E_value
            end
        end
    end
    return mean(energy_record), jackknife_energy_test(energy_record,n_loop,n_hist,sample_freq)
end

begin
    n_loops = 1000
    energy_list = []
    energy_error_list = []
    k_values = collect(range(0, 0.5, length=10))
    for k in k_values
        energy_vals, energy_error = dispersion_relation(1.0, k, n_loops)
        append!(energy_list, energy_vals)
        append!(energy_error_list, energy_error)
    end
end

begin
    quadratic(t, p) = p[1].+p[2].*t.+p[3].*t.*t
    p0=[0.0,-1.0,1.0]
    w=1 ./sqrt.(energy_error_list)
    fit = curve_fit(quadratic, k_values[1:10], energy_list[1:10], w,p0)
    println(1/(fit.param[3]*2))
    plot(k_values, energy_list, yerr = sqrt.(energy_error_list), label = "Dispersion")
    plot!(k_values, quadratic(k_values, fit.param), label = "Fit", xlabel="k", ylabel="E(k)")
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