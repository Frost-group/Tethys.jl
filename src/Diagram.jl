using Statistics
using Plots
using Distributions
using LinearAlgebra
using StaticArrays
using StructArrays

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

"""
    Arc(x...)

Type for storing the arc (phonon propagator) information, `x...`.
"""

mutable struct Arc <:Propagator
    q :: MVector{3,Float64}             # Phonon momentum
    period :: MVector{2,Float64}        # Arc length in the tau axis (scaled)
    ω :: Int64                          # Frequency
    index_in :: Int64                   # Index of the left vertex of the arc
    index_out :: Int64                  # Index of the right vertex of the arc
end


"""
    Line(x...)

Type for storing the line (electron propagator) information, `x...`.
"""

mutable struct Line <:Propagator
    k :: MVector{3,Float64}             # Electron momentum
    period :: MVector{2,Float64}        # Length of the line segment (scaled)
    index::Int64                        # Index of the line
    covered::Bool                       # True if the line is completely covered by a single arc
end


"""
    Diagram(x...)

Construtor for the complete diagram configuration.
"""

mutable struct Diagram

    mass :: Int64                           # Electron mass, set to 1
    μ :: Float64                            # Chemical potential                         
    ω :: Int64                              # Phonon frequency, set to 1
    α :: Float64                            # Coupling strength

    p :: MVector{3,Float64}                 # Total diagram momentum
    τ :: Float64                            # Total diagram length
    record_τ::Array{Float64,1}              # Array to record diagram lengths
    max_τ :: Float64                        # Maximum length of the diagram
    order :: Int64                          # Diagram order
    component :: Int64                      # Number of external phonons
    max_order :: Int64                      # Maximum diagram order
    
    arc_box :: Array{Arc,1}                 # Array of existing internal phonon arcs
    end_arc_box :: Array{Arc,1}             # Array of existing external phonon arcs
    line_box:: Array{Line,1}                # Array of electron line segments
    sign_box:: Array{Array{Int64, 1}, 1}    # Sign of the interaction vertices

    p_ins :: Float64                        # Probability of executing insert
    p_rem :: Float64                        # Probability of executing remove
    dispersion :: Float64                   # Value of the total dispersion of the diagram

    function Diagram(p::Real, max_τ::Real, max_order::Int64, mass::Int64, μ::Real, ω::Int64, α::Real, 
        p_ins::Float64=0.5, p_rem::Float64=0.5)
        τ = 20      # Diagram length set to 20
        new(mass, μ, ω, α, [p,0,0], τ, [], max_τ, 0, 0, max_order, [], [],
        [Line([p,0,0], [0, 1.0], 1,false)], [[0,0]],p_ins, p_rem,-τ*(norm(p)^2/(2*mass)-μ))
    end
end



"""
    green_zero(diagram)

Calculates the zeroth order diagram distribution

"""
function green_zero(diagram::Diagram)
    p=diagram.p
    τ=diagram.τ
    μ=diagram.μ
    m=diagram.mass

    return exp(-τ*(norm(p)^2/(2m)-μ))
end

function green_zero(line::Line, m::Int64, μ::Float64)
    p=line.k
    τ=line.period[2]-line.period[1]

    return exp(-τ*(norm(p)^2/(2m)-μ))
end


"""
    fast_norm(MVector)

Slightly faster way of calculating the norm of a vector

"""
function fast_norm(A::MVector{3,Float64})
    x = zero(eltype(A))
    for v in A
      x += v * v
    end
    x
  end


"""
  dispersion(line, m, μ)

Calculates the dispersion of a single line segment

"""
function dispersion(line::Line, m::Int64, μ::Float64)
    p=line.k
    τ=line.period[2]-line.period[1]

    return -τ*(norm(p)^2/(2m)-μ)
    #return -τ*(fast_norm(p)/(2m)-μ)
end

function p_dispersion(line::Line)
    p=line.k
    τ=line.period[2]-line.period[1]

    return p .*τ
end


"""
  dispersion(arc, ω)

Calculates the dispersion of an arc segment

"""
function dispersion(arc::Arc,ω::Int64)

    τ=arc.period[2]-arc.period[1]
    τ=line.period[2]-line.period[1]

    return -ω*(arc.period[2]-arc.period[1])
end


"""
  dispersion_vec(line_array, m, μ)

Vectorised version of the line dispersion calculation

"""
function dispersion_vec(line_array::StructArray{Line}, m::Int64, μ::Float64)
    p=line_array.k
    τ=getindex.(line_array.period,2)-getindex.(line_array.period,1)

    return -τ.*((p.⋅p)./(2m).-μ)
end

"""
  arc_dispersion(arc, ω)

Calculates the dispersion for an internal phonon propagator

"""
function arc_dispersion(arc::Arc, ω::Int64)
    return -ω*(arc.period[2]-arc.period[1])
end

"""
  arc_dispersion(arc, ω)

Calculates the dispersion for an external phonon propagator

"""
function end_arc_dispersion(arc::Arc, ω::Int64, τ::Float64)
    return -ω*(τ+arc.period[2]-arc.period[1])
end

"""
  phonon_propagator(arc)

Calculates the phonon propagator itself

"""
function phonon_propagator(arc::Arc)
    p=arc.q
    τ=arc.period[2]-arc.period[1]
    ω =arc.ω 

    return exp(-ω*τ)/norm(p)^2
end

function measure(diagram::Diagram)
    return [diagram.order,diagram.τ]
end

"""
  check_index(diagram)

Checks and displays the index of the line segments

"""
function check_index(diagram::Diagram)
    index_box=[]
    for index in 1:length(diagram.line_box)
        line=diagram.line_box[index]
        append!(index_box,line.index)
    end

    println("index_box",index_box)
end

"""
  check_arcindex(diagram)

Checks and displays the index of the arc segments

"""
function check_arcindex(diagram::Diagram)
    index_box=[]
    for index in 1:length(diagram.arc_box)
        arc=diagram.arc_box[index]
        push!(index_box,[arc.index_in,arc.index_out])
    end

    println("index_box",index_box)
end

"""
  check_timeorder(diagram)

Checks and displays the timeordering of the line segments

"""
function check_timeorder(diagram::Diagram)
    index_box=[]
    for index in 1:length(diagram.line_box)
        line=diagram.line_box[index]
        push!(index_box,[line.period[1],line.period[2]])
    end

    println("index_box",index_box)
    return index_box
end

"""
  check_timeorder(diagram)

Displays the electron momentum of the line segments

"""
function check_k(diagram::Diagram)
    index_box=[]
    for index in 1:length(diagram.line_box)
        line=diagram.line_box[index]
        push!(index_box,[line.k[1],line.k[2],line.k[3]])
    end

    println("index_box",index_box)
    return index_box
end
