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
        new(mass, μ, ω, α, [p_x,0,0], p_grid, τ, max_τ, 0, max_order, [], [Line([p_x,0,0], [0, τ], mass, μ)], p_ins, p_rem)
    end

end


abstract type Propagator end


mutable struct Arc <:Propagator
    q :: Array{Float64}
    period :: Array{Float64}
    ω :: Float64
    k_in :: Array{Float64}
    k_out :: Array{Float64}
end


mutable struct Line <:Propagator
    k :: Array{Float64}
    period :: Array{Float64}
    mass::Float64
    μ::Float64
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
        return
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
    new_arc=Arc(q,[τ_1,τ_2],ω,line.k,line.k)

    time=[τ_L,τ_1,τ_2,τ_R]
    k_box=[line.k,q,line.k]
    line_to_add=[]

    w_x=green_zero(line)
    w_y=phono_propagator(new_arc)*α_squared/(2*pi)^3

    for i in 1:3
        new_line=Line(k_box[i] ,[time[i],time[i+1]], m, μ)
        push!(line_to_add,new_line)
        w_y*=green_zero(new_line)
    end

    p_x_y=diagram.p_ins/(τ_R-τ_L)*ω*exp(-ω*arc_T)
    p_x_y*=exp(-norm(q)^2/(2m)*arc_T)/(2pi*m/arc_T)^1.5

    p_y_x=diagram.p_rem/(length(diagram.arc_box)+1)
    
    r=w_y*p_y_x/(w_x*p_x_y)
    println("r=",r)

    if r<rand()
        return false
    else
        diagram.order+=1
        push!(diagram.arc_box,new_arc)
        deleteat!(line_box, index)

        for i in 1:3
            insert!(line_box, index, line_to_add[4-i])
        end

        return true
    end
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
    loop=100
    for i in 1:loop
        insert_arc!(diagram_a)
        println("order is ",diagram_a.order)
    end
end

