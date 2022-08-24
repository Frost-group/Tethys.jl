using Statistics
using Plots
using Distributions
using LinearAlgebra

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

    p_ins :: Float64
    p_rem :: Float64

    function Diagram(p_max::Real, max_τ::Real, max_order::Int64, mass::Real, μ::Real, ω::Real, α::Real, 
        p_ins::Float64=0.5, p_rem::Float64=0.5, grid_num::Int64=1000)
        τ = rand()*max_τ
        p_grid=collect(LinRange(0, p_max, grid_num))
        p_x=rand(p_grid)
        new(mass, μ, ω, α, [p_x,0,0], p_grid, τ, max_τ, 0, max_order, [], [Line([p_x,0,0], [0, τ], mass, μ, 1)], p_ins, p_rem)
    end

end


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

function phono_propagator(arc::Arc)
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
    line_to_add=[]

    w_x=green_zero(line)
    w_y=phono_propagator(new_arc)*α_squared/(2*pi)^3

    for i in 1:3
        new_line=Line(k_box[i] ,[time[i],time[i+1]], m, μ, index+i-1)
        push!(line_to_add,new_line)
        w_y*=green_zero(new_line)
    end

    p_x_y=diagram.p_ins/(τ_R-τ_L)*ω*exp(-ω*arc_T)
    p_x_y*=exp(-norm(q)^2/(2m)*arc_T)/(2pi*m/arc_T)^1.5

    p_y_x=diagram.p_rem/(order+1)
    
    r=w_y*p_y_x/(w_x*p_x_y)
    # println("insert_r=",r)

    if r<rand()
        return false
    else
        diagram.order+=1
        deleteat!(line_box, index)

        for i in 1:3
            insert!(line_box, index, line_to_add[4-i])
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
    new_line=Line(line_in.k ,[τ_L,τ_R], m, μ, index_in)

    w_x=green_zero(new_line)
    w_y=phono_propagator(arc)*α_squared/(2*pi)^3

    for i in 1:3
        w_y*=green_zero(line_to_rem[i])
    end

    τ_1=copy(arc.period[1])
    τ_2=copy(arc.period[2])
    arc_T=τ_2-τ_1
    p_x_y=diagram.p_ins/(τ_R-τ_L)*ω*exp(-ω*arc_T)
    p_x_y*=exp(-norm(q)^2/(2m)*arc_T)/(2pi*m/arc_T)^1.5

    p_y_x=diagram.p_rem/order

    r=(w_x*p_x_y)/(w_y*p_y_x)
    # println("remove_r=",r)

    if r<rand()
        return false
    else
        diagram.order-=1
        
        deleteat!(arc_box, index)
        deleteat!(line_box, index_in:index_out)
        insert!(line_box, index_in, new_line)

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
            # elseif arc.index_in>index_in && arc.index_out<index_out
            #     arc.index_in-=1
            #     arc.index_out-=1
            end
        end

        return true
    end
end

function swap_arc!(diagram::Diagram)

    order=diagram.order
    m=diagram.mass
    μ=diagram.μ
    ω=diagram.ω
    α=diagram.α
    α_squared=2pi*α*sqrt(2)

    if order<2
        return false
    end

    arc_box=diagram.arc_box
    sort!(arc_box, by = x -> x.index_in)
    index=rand(1:length(arc_box)-1)
    arc_1=arc_box[index]
    arc_2=arc_box[index+1]

    if (arc_1.index_out-arc_1.index_in)>2 || (arc_2.index_out-arc_2.index_in)>2 #|| arc_1.index_out != arc_2.index_in
        return false
    end

    τ_a=arc_1.period[1]
    τ_1=arc_1.period[2]
    τ_2=arc_2.period[1]
    τ_b=arc_2.period[2]
    q1=arc_1.q
    q2=arc_2.q
    line_box=diagram.line_box
    line=line_box[arc_1.index_out]

    w_x=green_zero(line)*phono_propagator(arc_1)*phono_propagator(arc_2)

    new_arc1=Arc(q1,[τ_a,τ_2],ω,arc_1.index_in,arc_1.index_out+1)
    new_arc2=Arc(q2,[τ_1,τ_b],ω,arc_2.index_in-1,arc_2.index_out)
    new_line=Line(line.k-q1-q2 ,[τ_1,τ_2], m, μ, line.index)

    w_y=green_zero(new_line)*phono_propagator(new_arc1)*phono_propagator(new_arc2)

    r=w_y/w_x

    if r<rand()
        return false
    else
        deleteat!(line_box, line.index)
        insert!(line_box, line.index, new_line)
        deleteat!(arc_box, index:index+1)
        insert!(arc_box, index, new_arc1)
        insert!(arc_box, index+1, new_arc2)
        return true
    end
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

begin
    p_max=10; max_τ=20; max_order=50; mass=1; μ=1.6; ω=1; α=1.2
    diagram_a=Diagram(p_max, max_τ, max_order, mass, μ, ω, α)

    loop=10
    record=[]
    for i in 1:loop
        τ_update!(diagram_a)
        append!(record,diagram_a.τ)
    end
    histogram(record)
end

begin
    println("begin")
    diagram_a=Diagram(p_max, max_τ, max_order, mass, μ, ω, α)
    loop=100
    for i in 1:loop

        if insert_arc!(diagram_a)
            println("insert")
            check_index(diagram_a)
            # check_arcindex(diagram_a)
            # check_timeorder(diagram_a)
        end
        if remove_arc!(diagram_a)
            println("remove")
            # check_arcindex(diagram_a)
            check_index(diagram_a)
            # check_timeorder(diagram_a)
        end
        if swap_arc!(diagram_a)
            println("swap")
            # check_arcindex(diagram_a)
            check_index(diagram_a)
            # check_timeorder(diagram_a)
        end

        # if insert_arc!(diagram_a)
        #     check_arcindex(diagram_a)
        #     check_timeorder(diagram_a)
        #     println("order is ",diagram_a.order)
        # end
        # println("order is ",diagram_a.order)
    end
end

