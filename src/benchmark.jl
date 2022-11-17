# include("Diagram.jl")
# include("update.jl")
# include("measure.jl")
# include("Tethys.jl")
import Tethys
using Random
using Logging
using BenchmarkTools

function insert_arc_benchmark(α, μ, p=0)
    max_τ=30; max_order=500; mass=1; ω=1;
    diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)
    order=diagram.order
    m=diagram.mass
    μ=diagram.μ
    ω=diagram.ω
    α=diagram.α
    α_squared=2pi*α*sqrt(2)
    insert_arc!(diagram,order,m,μ,ω,α_squared)
    order=diagram.order
    insert_arc!(diagram,order,m,μ,ω,α_squared)
    order=diagram.order
    insert_arc!(diagram,order,m,μ,ω,α_squared)
end

function insert_arc_2_benchmark(α, μ, p=0)
    max_τ=30; max_order=500; mass=1; ω=1;
    regime=Diff_more()
    diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)
    insert_arc_2!(diagram)
    insert_arc_2!(diagram)
    insert_arc_2!(diagram)
end

begin
    Random.seed!(1234)
    t = @benchmark insert_arc_benchmark(5.0, -6)
end

begin
    Random.seed!(1234)
    t = @benchmark insert_arc_2_benchmark(5.0, -6)
end

function testing_arr(arr)
    q,w,e,r = arr

end
begin
    t = @benchmark testing_arr([1,2,3,4])
end

begin
    n_loop=10
    num_mea=1; regime=Diff_more();
    p=0; max_τ=30; max_order=500; mass=1; ω=1;


    linear(t, p) = p[1].-p[2].*t
    bin_width=max_τ/300
    min_time=Int(div(5,bin_width,RoundUp))
    max_time=Int(div(12,bin_width,RoundUp))

    α=5
    μ=-6
    directory = "F://data"
    store_data = false
    hist=Hist_Record(300,max_τ,max_order)
    diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)
    diagram,hist,green_record,zero_record,green_func,variance=hist_measure!(diagram,hist,directory,n_loop,store_data)
    
    @info "Diagram Order: "*string(diagram.order)
    t = @benchmark insert_arc!(diagram,order,mass,μ,ω,α_squared)
end

function insert_arc_2!(diagram::Diagram)

    order=diagram.order
    m=diagram.mass
    μ=diagram.μ
    ω=diagram.ω
    α=diagram.α
    α_squared=2pi*α*sqrt(2)
    τ=diagram.τ

    if order+1>diagram.max_order
        return false
    end

    line_box=diagram.line_box
    index_in=rand(1:length(line_box))
    line=line_box[index_in]
    τ_L=deepcopy(line.period[1])
    τ_R=deepcopy(line.period[2])
    k_in=deepcopy(line.k)

    τ_1=rand(Uniform(τ_L,τ_R))
    τ_2=τ_1-log(rand())/ω 

    if τ_2 > τ
        return false
    end

    line.period[1]=τ_1
    arc_T=τ_2-τ_1
    q=rand(Normal(0,sqrt(m/arc_T)),3)
    
    w_x=1
    w_y=1
    total_dis=0
    τ_R_2=0
    index_out=0
    k_out=0

    #not set covered yet
    for i in index_in:2order+1
        line_tem=line_box[i]
        if line_tem.period[2]<τ_2
            total_dis+=dispersion(line_tem)
            #w_x*=green_zero(line_tem)
            line_tem.k-=q
            line_tem.index+=1
            total_dis-=dispersion(line_tem)
            #w_x/=green_zero(line_tem)
            continue
        else
            k_out=deepcopy(line_tem.k)
            τ_R_2=deepcopy(line_tem.period[2])
            line_tem.period[2]=τ_2
            total_dis+=dispersion(line_tem)
            #w_x*=green_zero(line_tem)
            line_tem.k-=q
            line_tem.index+=1
            total_dis-=dispersion(line_tem)
            #w_x/=green_zero(line_tem)
            index_out=i+2
            break
        end
    end

    # println("index_in:",index_in)
    # println("index_out",index_out-2)
    new_arc=Arc(q,[τ_1,τ_2],ω,index_in,index_out)

    #w_y*=phonon_propagator(new_arc)*α_squared/(2*pi)^3
    p_x_y=diagram.p_ins/(2order+1)/(τ_R-τ_L)#*ω*exp(-ω*arc_T)#/(1-exp(-ω*(τ-τ_1)))#
    p_x_y*=exp(-norm(q)^2/(2m)*arc_T)/(2pi*m/arc_T)^1.5#/normalization(τ,τ_R,τ_L,ω)#*weighting/exp(-ω*τ_1)
    # println("ratio=",line_length/(τ-τ_1))
    p_y_x=diagram.p_rem/(order+1)
    #r=w_y*p_y_x/(w_x*p_x_y)
    r=α_squared*p_y_x/(exp(total_dis)*p_x_y*(2*pi)^3*norm(q)^2)
    # println("normal",normalization(τ,τ_R,τ_L,ω))
    # println("insert_r=",r)
    

    if r<rand()
        line_box[index_in].period[1]=τ_L
        line_box[index_out-2].period[2]=τ_R_2
        for i in index_in:index_out-2
            line_tem=line_box[i]
            line_tem.index=i
            line_tem.k+=q
        end

        return false
    else
        diagram.order+=1

        if index_out-index_in==2
            line_box[index_in].covered=true
        else
            line_box[index_in].covered=false
            line_box[index_out-2].covered=false
        end

        line_tem=Line(k_in,[τ_L,τ_1], m, μ, index_in, false)
        insert!(line_box, index_in, line_tem)
        line_tem=Line(k_out,[τ_2,τ_R_2], m, μ, index_out, false)
        insert!(line_box, index_out, line_tem)

        sign_box=diagram.sign_box
        if index_out-index_in==2
            sign_to_add=[[copy(sign_box[index_in][1]),-1],[-1,1],[1,copy(sign_box[index_in][2])]]
            deleteat!(sign_box, index_in)
            for i in 1:3
                insert!(sign_box, index_in, sign_to_add[4-i])
            end
            
        else

            sign_to_add=[[copy(sign_box[index_out-2][1]),1],[1,copy(sign_box[index_out-2][2])]]
            deleteat!(sign_box, index_out-2)
            for i in 1:2
                insert!(sign_box, index_out-2, sign_to_add[3-i])
            end

            sign_to_add=[[copy(sign_box[index_in][1]),-1],[-1,copy(sign_box[index_in][2])]]
            deleteat!(sign_box, index_in)
            for i in 1:2
                insert!(sign_box, index_in, sign_to_add[3-i])
            end
        end

        if length(line_box)>=index_out+1
            for i in index_out+1:length(line_box)
                line_box[i].index=i
            end
        end

        for arc in diagram.arc_box
            if arc.index_out<=index_in
                continue
            elseif arc.index_in>=index_out-2
                arc.index_in+=2
                arc.index_out+=2
            elseif arc.index_in<index_in && arc.index_out>index_out-2
                arc.index_out+=2
            elseif arc.index_in>=index_in && arc.index_out<=index_out-2
                arc.index_in+=1
                arc.index_out+=1
            elseif index_in<=arc.index_in<index_out-2
                arc.index_in+=1
                arc.index_out+=2
            elseif index_in<arc.index_out<=index_out-2
                arc.index_out+=1
            end
        end

        push!(diagram.arc_box,new_arc)

        return true
    end
end