include("Diagram.jl")
include("update.jl")
include("zero_update.jl")
include("measure.jl")
using Random
using LsqFit
using JLD2
using FFTW
using Logging
using Statistics
using Dates
using CSV, DataFrames
using JLD2,FileIO
using PolaronMobility

mutable struct Estimators_Record

    max_τ::Float64
    max_order::Int64
    sample_freq::Int64
    n_loop::Int64
    n_hist::Int64
    weight_box::Array{Float64,1}
    order_box::Array{Float64,1}
    energy_record::Array{Float64,1}
    energy_mean::Array{Float64,1}
    p_record::Array{Float64,1}
    mass_mean::Array{Float64,1}

    function Estimators_Record(max_τ::Float64,max_order::Int64,sample_freq::Int64,n_loop::Int64, n_hist::Int64)
        weight_box = zeros(max_order)
        order_box = zeros(max_order)
        energy_record=zeros(Int64(n_loop*n_hist/sample_freq))
        energy_mean=zeros(Int64(n_loop*n_hist/sample_freq))
        p_record=zeros(Int64(n_loop*n_hist/sample_freq))
        mass_mean=zeros(Int64(n_loop*n_hist/sample_freq))
        new(max_τ,max_order,sample_freq,n_loop,n_hist,weight_box,
        order_box,energy_record,energy_mean,p_record,mass_mean)
    end
end

"""
    VMC_energy(α)
Returns the zero-temperature ground-state energy of the polaron for a material with multiple phonon branches.
From the PolaronMobility Julia package which uses the Feynman's variational method. 
See: https://github.com/jarvist/PolaronMobility.jl
"""
function VMC_energy(α)
    v,w = feynmanvw(α)
    return F(v, w, α)
end

function polaron_effective_mass(α)
    v,w = feynmanvw(α)
    massF(τ,v,w)=(abs(w^2 * τ + (v^2-w^2)/v*(1-exp(-v*τ))))^-1.5 * exp(-τ) * τ^2
    intF(v,w,α)=(1/3)*π^(-0.5) * α*v^(3) * quadgk(τ->massF(τ,v,w),0,Inf)[1]

    return 1+intF(v,w,α)
end

function jackknife_energy(energy_samples)
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

function loop_verbose(loop_number, frequency=10)
    if loop_number%frequency != 0
        return false
    else
        @info "$(Dates.format(now(), "yyyy-mm-dd HH:MM:SS")) loop.number:"*string(loop_number)
    end
end

function initialise_diagram(α::Float64, p::Float64, max_τ::Float64, max_order::Int64, μ::Float64=0.0)
    ω = 1
    mass = 1
    return Diagram(p, max_τ, max_order, mass, μ, ω, α)
end

function simulate!(diagram::Diagram,estimators::Estimators_Record, swap_arc=false, store_data=false,
    p_ins=0.2,p_rem=0.2,p_from_0=1)

    real_normalized=[p_ins,p_rem]
    real_normalized/=sum(real_normalized)
    fake_normalized=[p_from_0]
    fake_normalized/=sum(fake_normalized)
    real_cumsum=cumsum(real_normalized)
    diagram.p_ins=real_normalized[1]
    diagram.p_rem=real_normalized[2]
    weight_box = estimators.weight_box
    order_box= estimators.order_box
    sample_freq = estimators.sample_freq
    energy_record=estimators.energy_record
    energy_mean=estimators.energy_mean
    p_record=estimators.p_record
    mass_mean=estimators.mass_mean

    n_loop=estimators.n_loop
    n_hist=estimators.n_hist
    dia_order=diagram.order
    m=diagram.mass
    μ=diagram.μ
    ω=diagram.ω
    α=diagram.α
    α_squared=2pi*α*sqrt(2)
    τ=diagram.τ

    println("begin")
    for j in 1:n_loop
        loop_verbose(j)
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

            dia_order=diagram.order

            if !result && swap_arc
                result = swap_arc!(diagram)
                #result = resample_arc!(diagram,dia_order,m,μ,ω,α_squared)
                #swap_arc!(diagram)
            end

            
            #scale!(diagram,dia_order,m,μ,ω,1)
            #update_arcp!(diagram,dia_order,m,μ)

            component=diagram.component
            weight_box[component+1]+=1
            order_box[dia_order+1]+=1

            #while argmax(weight_box) <= 40
            #    println(argmax(weight_box))
            #    insert_arc!(diagram,dia_order,m,μ,ω,α_squared)
            #    zero_remove_arc!(diagram,dia_order,m,μ,ω,α_squared)
            #    dia_order=diagram.order
            #    component=diagram.component
            #    weight_box[component+1]+=1
            #    order_box[dia_order+1]+=1
            #end

            if mod(i,sample_freq) == 0
                estimator_index = Int64(((j-1)*n_hist+i)/sample_freq)
                E_value=energy(diagram)
                p_value=mass_estimator(diagram)
                energy_record[estimator_index] = E_value
                p_record[estimator_index] = p_value
                if estimator_index == 1  
                    energy_mean[estimator_index] = mean(energy_record[1:estimator_index])
                    #mass_mean[estimator_index] = 1/(1-mean(p_record[1:estimator_index])*τ/3)
                    #mass_mean[estimator_index] = result[2]
                    mass_mean[estimator_index] = τ/(2*dia_order+1)
                else
                    energy_mean[estimator_index] = (E_value + energy_mean[estimator_index-1]*(estimator_index-1))/estimator_index
                    p_mean = (p_value + (3/τ)*(1-1/mass_mean[estimator_index-1])*(estimator_index-1))/estimator_index
                    #mass_mean[estimator_index] = 1/(1-p_mean*τ/3)
                    mass_mean[estimator_index] = τ/(2*dia_order+1)
                    #mass_mean[estimator_index] = diagram.line_box[1].period[2]-diagram.line_box[1].period[1]
                    #mass_mean[estimator_index] = result[2]
                end
                
            end
        end
    end
    return diagram, estimators, jackknife_energy(energy_record)
end