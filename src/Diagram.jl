using Statistics
using Plots
using Distributions
using LinearAlgebra

abstract type Propagator end


mutable struct Arc <:Propagator
    q :: Array{Float64}
    period :: Array{Float64}
    ω :: Float64
    index_in :: Int64
    index_out :: Int64
end


mutable struct Line <:Propagator
    k :: Array{Float64}
    period :: Array{Float64}
    mass::Float64
    μ::Float64
    index::Int64
    covered::Bool
end

mutable struct Diagram

    mass:: Float64
    μ::Float64
    ω::Float64
    α::Float64

    p :: Array{Float64}
    p_grid :: Array{Float64}
    τ :: Float64
    max_τ :: Float64
    order :: Int64
    max_order :: Int64
    
    arc_box :: Array{Arc}
    line_box:: Array{Line}
    sign_box:: Array{Array{Int64, 1}, 1}

    p_ins :: Float64
    p_rem :: Float64

    function Diagram(p_max::Real, max_τ::Real, max_order::Int64, mass::Real, μ::Real, ω::Real, α::Real, 
        p_ins::Float64=0.5, p_rem::Float64=0.5, grid_num::Int64=1000)
        τ = rand()*max_τ
        p_grid=collect(LinRange(0, p_max, grid_num))
        p_x=rand(p_grid)
        new(mass, μ, ω, α, [p_x,0,0], p_grid, τ, max_τ, 0, max_order, [], 
        [Line([p_x,0,0], [0, τ], mass, μ, 1,false)], [[0,0]],p_ins, p_rem)
    end

end

function green_zero(diagram::Diagram)
    p=diagram.p
    τ=diagram.τ
    μ=diagram.μ
    m=diagram.mass

    return exp(-τ*(norm(p)^2/(2m)-μ))
end

function green_zero(line::Line)
    p=line.k
    τ=line.period[2]-line.period[1]
    μ=line.μ
    m=line.mass

    return exp(-τ*(norm(p)^2/(2m)-μ))
end

function phonon_propagator(arc::Arc)
    p=arc.q
    τ=arc.period[2]-arc.period[1]
    ω =arc.ω 

    return exp(-ω*τ)/norm(p)^2
end

function p_update!(diagram::Diagram)

    if diagram.order != 0
        return 
    end

    p_sample=diagram.p_grid
    p_new=[rand(p_sample)*rand([-1,1]),0,0]
    p_old=copy(diagram.p)
    old_green=green_zero(diagram)
    diagram.p=p_new
    new_green=green_zero(diagram)

    if new_green/old_green >= rand()
        return true
    else
        diagram.p=p_old
        return false
    end
end

function τ_update!(diagram::Diagram)

    if diagram.order != 0
        return 
    end

    p=diagram.p
    μ=diagram.μ
    m=diagram.mass
    dispersion=norm(p)^2/(2m)-μ

    τ_new=-log(rand())/abs(dispersion)

    if τ_new>=diagram.max_τ
        return false
    else
        diagram.τ=τ_new
        return true
    end
end

function insert_arc!(diagram::Diagram)

    order=diagram.order
    m=diagram.mass
    μ=diagram.μ
    ω=diagram.ω
    α=diagram.α
    α_squared=2pi*α*sqrt(2)

    if order+1>diagram.max_order
        return false
    end

    line_box=diagram.line_box
    index=rand(1:length(line_box))
    line=line_box[index]
    τ_L=line.period[1]
    τ_R=line.period[2]

    τ_1=rand(Uniform(τ_L,τ_R))
    τ_2=τ_1-log(rand())/diagram.ω 

    if τ_2 > τ_R
        #println("no correct τ_2")
        return false
    end

    arc_T=τ_2-τ_1

    q=rand(Normal(0,sqrt(m/arc_T)),3)
    new_arc=Arc(q,[τ_1,τ_2],ω,index,index+2)

    time=[τ_L,τ_1,τ_2,τ_R]
    k_box=[line.k,line.k-q,line.k]
    covered=[false,true,false]
    line_to_add=[]

    w_x=green_zero(line)
    w_y=phonon_propagator(new_arc)*α_squared/(2*pi)^3

    for i in 1:3
        new_line=Line(k_box[i] ,[time[i],time[i+1]], m, μ, index+i-1, covered[i])
        push!(line_to_add,new_line)
        w_y*=green_zero(new_line)
    end

    p_x_y=diagram.p_ins/(τ_R-τ_L)*ω*exp(-ω*arc_T)
    p_x_y*=exp(-norm(q)^2/(2m)*arc_T)/(2pi*m/arc_T)^1.5

    p_y_x=diagram.p_rem/(order+1)
    
    r=w_y*p_y_x/(w_x*p_x_y)#*(2order+1)
    # println("insert_r=",r)

    if r<rand()
        return false
    else
        diagram.order+=1
        deleteat!(line_box, index)
        sign_box=diagram.sign_box
        sign_to_add=[[copy(sign_box[index][1]),-1],[-1,1],[1,copy(sign_box[index][2])]]
        deleteat!(sign_box, index)

        for i in 1:3
            insert!(line_box, index, line_to_add[4-i])
            insert!(sign_box, index, sign_to_add[4-i])
        end
        # println("insert index is ",[index,index+2])
        if length(line_box)>=index+3
            for i in index+3:length(line_box)
                line_box[i].index=i
            end
        end

        for arc in diagram.arc_box
            if arc.index_out<=index
                continue
            elseif arc.index_in>=index
                arc.index_in+=2
                arc.index_out+=2
            elseif arc.index_in<index && arc.index_out>index
                arc.index_out+=2
            end
        end

        push!(diagram.arc_box,new_arc)

        return true
    end
end

function remove_arc!(diagram::Diagram)

    order=diagram.order
    m=diagram.mass
    μ=diagram.μ
    ω=diagram.ω
    α=diagram.α
    α_squared=2pi*α*sqrt(2)

    if order-1<1
        return false
    end
    
    arc_box=diagram.arc_box
    index=rand(1:length(arc_box))
    arc=arc_box[index]
    index_in=arc.index_in
    index_out=arc.index_out
    q=arc.q

    if index_out-index_in>2
        return false
    end

    line_box=diagram.line_box
    line_in=line_box[index_in]
    line_out=line_box[index_out]
    line_to_rem=[line_box[index_in],line_box[index_in+1],line_box[index_in+2]]

    τ_L=line_in.period[1]
    τ_R=line_out.period[2]
    new_line=Line(line_in.k ,[τ_L,τ_R], m, μ, index_in,false)

    w_x=green_zero(new_line)
    w_y=phonon_propagator(arc)*α_squared/(2*pi)^3

    for i in 1:3
        w_y*=green_zero(line_to_rem[i])
    end

    τ_1=copy(arc.period[1])
    τ_2=copy(arc.period[2])
    arc_T=τ_2-τ_1
    p_x_y=diagram.p_ins/(τ_R-τ_L)*ω*exp(-ω*arc_T)
    p_x_y*=exp(-norm(q)^2/(2m)*arc_T)/(2pi*m/arc_T)^1.5

    p_y_x=diagram.p_rem/order

    r=(w_x*p_x_y)/(w_y*p_y_x)#/(2order-1)
    # println("remove_r=",r)

    if r<rand()
        return false
    else
        diagram.order-=1
        
        deleteat!(arc_box, index)
        deleteat!(line_box, index_in:index_out)
        insert!(line_box, index_in, new_line)
        covered=false

        sign_box=diagram.sign_box
        sign_to_add=[sign_box[index_in][1],sign_box[index_out][2]]
        deleteat!(sign_box, index_in:index_out)
        insert!(sign_box, index_in, sign_to_add)
        

        if length(line_box)>=index_in+1
            for i in index_in+1:length(line_box)
                line_box[i].index=i
            end
        end

        for arc in diagram.arc_box
            if arc.index_out<=index_in
                continue
            elseif arc.index_in>=(index_in+2)
                arc.index_in-=2
                arc.index_out-=2
            elseif arc.index_in<index_in && arc.index_out>index_in
                arc.index_out-=2
                covered=true
            # elseif arc.index_in>index_in && arc.index_out<index_out
            #     arc.index_in-=1
            #     arc.index_out-=1
            end
        end

        if covered
            line_box[index_in].covered=covered
        end

        return true
    end
end

function arc_judge(arc::Arc,sign::Int64,bound::Bool,index::Int64)
    # bound true is right, false is left
    if bound 
        if sign == 1
            return arc.index_out-1 == index
        else
            return arc.index_in == index
        end
    else
        if sign == 1
            return arc.index_out == index
        else
            return arc.index_in+1 == index
        end
    end
end

function swap_arc!(diagram::Diagram)

    order=diagram.order
    m=diagram.mass
    μ=diagram.μ
    ω=diagram.ω

    if order<2
        return false
    end

    line_box=diagram.line_box
    line_index=rand(2:2*order)
    chosen_line=line_box[line_index]

    if chosen_line.covered
        return false
    end

    sign=diagram.sign_box[line_index]
    arc_box=diagram.arc_box
    left_check=false
    right_check=false
    left_index=0
    right_index=0

    for i in 1:order
        arc=arc_box[i]
        if !left_check
            if arc_judge(arc,sign[1],false,line_index)
                left_index=i
                left_check=true
            end
        end

        if !right_check
            if arc_judge(arc,sign[2],true,line_index)
                right_index=i
                right_check=true
            end
        end

        if right_check && left_check
            break
        end
    end

    arc_l=arc_box[left_index]
    arc_r=arc_box[right_index]

    q1=arc_l.q
    q2=arc_r.q

    if sign[1] == 1
        new_arc_l=Arc(q1,[arc_l.period[1],chosen_line.period[2]],ω,arc_l.index_in,line_index+1)
    else
        new_arc_l=Arc(q1,[chosen_line.period[2],arc_l.period[2]],ω,line_index,arc_l.index_out)
    end

    if sign[2] == 1
        new_arc_r=Arc(q2,[arc_r.period[1],chosen_line.period[1]],ω,arc_r.index_in,line_index)
    else
        new_arc_r=Arc(q2,[chosen_line.period[1],arc_r.period[2]],ω,line_index-1,arc_r.index_out)
    end

    new_line=Line(chosen_line.k-sign[1]*q1+sign[2]*q2 ,chosen_line.period, m, μ, line_index,false)
    w_x=green_zero(chosen_line)*phonon_propagator(arc_l)*phonon_propagator(arc_r)
    w_y=green_zero(new_line)*phonon_propagator(new_arc_l)*phonon_propagator(new_arc_r)

    r=w_y/w_x

    if r<rand()
        return false
    else
        deleteat!(line_box, line_index)
        insert!(line_box, line_index, new_line)
        sign_box=diagram.sign_box
        deleteat!(sign_box, line_index)
        insert!(sign_box, line_index, [sign[2],sign[1]])
        sign_box[line_index-1]=[sign_box[line_index-1][1],sign[2]]
        sign_box[line_index+1]=[sign[1],sign_box[line_index+1][2]]

        deleteat!(arc_box, left_index)
        insert!(arc_box, left_index, new_arc_l)
        deleteat!(arc_box, right_index)
        insert!(arc_box, right_index, new_arc_r)

        if new_arc_l.index_out-new_arc_l.index_in == 2
            line_box[line_index-1].covered=true
        end

        if new_arc_r.index_out-new_arc_r.index_in == 2
            line_box[line_index+1].covered=true

        end
        return true
    end
end

function extend!(diagram::Diagram)

    p=diagram.p
    μ=diagram.μ
    m=diagram.mass
    dispersion=norm(p)^2/(2m)-μ
    line_box=diagram.line_box
    order=diagram.order
    line_end=line_box[2*order+1]

    τ_new=line_end.period[1]-log(rand())/abs(dispersion)

    if τ_new>diagram.max_τ
        return false
    end

    line_end.period[2]=τ_new
    line_box[2*order+1]=line_end
    diagram.τ=τ_new

    return true

end

function measure(diagram::Diagram)
    return [diagram.order,diagram.τ]
end

function check_index(diagram::Diagram)
    index_box=[]
    for index in 1:length(diagram.line_box)
        line=diagram.line_box[index]
        append!(index_box,line.index)
    end

    println("index_box",index_box)
end

function check_arcindex(diagram::Diagram)
    index_box=[]
    for index in 1:length(diagram.arc_box)
        arc=diagram.arc_box[index]
        push!(index_box,[arc.index_in,arc.index_out])
    end

    println("index_box",index_box)
end

function check_timeorder(diagram::Diagram)
    index_box=[]
    for index in 1:length(diagram.line_box)
        line=diagram.line_box[index]
        push!(index_box,[line.period[1],line.period[2]])
    end

    println("index_box",index_box)
end

