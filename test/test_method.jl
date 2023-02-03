using Tethys
using Random
using LsqFit
using JLD2

begin
    n_loop=200
    num_samples=20
    n_hist=100000
    α=5
    μ=-5.6
    num_mea=1; regime=Diff_more();
    p=0; max_τ=30; max_order=500; mass=1; ω=1;
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
            extend!(diagram)
            dia_order=diagram.order
        end
    end
end