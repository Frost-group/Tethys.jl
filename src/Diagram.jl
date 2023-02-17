using Statistics
using Plots
using Distributions
using LinearAlgebra
using StaticArrays

abstract type Propagator end

abstract type Regime end

struct Diff_2 <: Regime
    function Diff_2()
        new()
    end
end

struct Diff_more <: Regime
    function Diff_more()
        new()
    end
end

mutable struct Arc <:Propagator
    q :: MVector{3,Float64}
    period :: MVector{2,Float64}
    ω :: Float64
    index_in :: Int64
    index_out :: Int64
end


mutable struct Line <:Propagator
    k :: MVector{3,Float64}
    period :: MVector{2,Float64}
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

    p :: MVector{3,Float64}
    τ :: Float64
    record_τ::Float64
    max_τ :: Float64
    order :: Int64
    max_order :: Int64
    
    arc_box :: Array{Arc}
    end_arc_box :: Array{Arc}
    line_box:: Array{Line}
    sign_box:: Array{Array{Int64, 1}, 1}

    p_ins :: Float64
    p_rem :: Float64

    function Diagram(p::Real, max_τ::Real, max_order::Int64, mass::Real, μ::Real, ω::Real, α::Real, 
        p_ins::Float64=0.5, p_rem::Float64=0.5)
        τ = 7.5#5.5/ω+rand()*max_τ
        new(mass, μ, ω, α, [p,0,0], τ,τ, max_τ, 0, max_order, [], [],
        [Line([p,0,0], [0, τ], mass, μ, 1,false)], [[0,0]],p_ins, p_rem)
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

function fast_norm(A::MVector{3,Float64})
    x = zero(eltype(A))
    for v in A
      x += v * v
    end
    x
  end

function dispersion(line::Line)
    p=line.k
    τ=line.period[2]-line.period[1]
    μ=line.μ
    m=line.mass

    return -τ*(norm(p)^2/(2m)-μ)
end

function phonon_propagator(arc::Arc)
    p=arc.q
    τ=arc.period[2]-arc.period[1]
    ω =arc.ω 

    return exp(-ω*τ)/norm(p)^2
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
