#using .Tethys
include("../src/Diagram.jl")
include("../src/update.jl")
include("../src/measure.jl")
using Random
using LsqFit
using JLD2

begin
    n_loop=1000
    num_samples=20
    n_hist=100000
    α=17
    μ=-34
    num_mea=1; regime=Diff_more();
    p=0; max_τ=300; max_order=2000; mass=1; ω=1;
end

begin
    Random.seed!(13653)
    diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)

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
                insert_arc!(diagram,dia_order,m,μ,ω,α_squared)
                diagram.p_ins=real_normalized[1]
            elseif  dia_order == 1
                if q<real_cumsum[1]
                    insert_arc!(diagram,dia_order,m,μ,ω,α_squared)
                else
                    diagram.p_ins=fake_normalized[1]
                    remove_arc!(diagram,dia_order,m,μ,ω,α_squared)
                    diagram.p_ins=real_normalized[1]
                end
            else
                if q<real_cumsum[1]
                    insert_arc!(diagram,dia_order,m,μ,ω,α_squared)
                else
                    remove_arc!(diagram,dia_order,m,μ,ω,α_squared)
                end
            end
            dia_order=diagram.order
            swap_arc!(diagram)
            extend!(diagram)
            τ = diagram.τ
            unnormalized_data[Int(div(τ, bin_width, RoundUp))]+=1
            if τ>6
                a = length(diagram.end_arc_box)
                weight_box[a+1]+=1
            end
        end
    end
end

begin
    plot(weight_box/sum(weight_box), xlims=(20,80))
end

begin
    time_points=collect(1:num_bins)*bin_width.-(bin_width/2)
    display(plot(time_points, log.(unnormalized_data), xlims=(0,50)))
    linear(t, p) = p[1].-p[2].*t
    min_time=Int(div(9,bin_width,RoundUp))
    max_time=Int(div(14,bin_width,RoundUp))
    time_points=time_points[min_time:max_time]
    p0=[0,(-α-1.26*(α/10)^2-μ)]
    y=log.(unnormalized_data[min_time:max_time])
    fit = curve_fit(linear, time_points, y, p0)
    println("energy:",fit.param[2]+μ)
    println("perturb:",-α-1.26*(α/10)^2)

end