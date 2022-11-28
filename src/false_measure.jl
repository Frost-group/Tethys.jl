include("Diagram.jl")
include("update.jl")
include("measure.jl")
using Random
using CSV, DataFrames
using JLD2


mutable struct Discrete_Record

    min_τ:: Float64
    max_τ::Float64
    time_points::Array{Float64}
    max_order::Int64
    unnormalized_data::Array{Int64,2}
    normalized_data::Array{Float64,2}

    function Discrete_Record(num_points::Int64,min_τ::Real,max_τ::Real,max_order::Int64)
        bin_width=(max_τ-min_τ)/(num_points-1)
        unnormalized_data=zeros(max_order+1, num_points)
        time_points=collect(0:num_points-1)*bin_width.+ min_τ
        new(min_τ, max_τ, time_points, max_order,unnormalized_data,unnormalized_data)
    end

end

function normalization(data,time_points,index,diagram::Diagram)
    p=diagram.p
    μ=diagram.μ
    m=diagram.mass
    zeroth_order=deepcopy(data[1,index])
    τ=time_points[index]
    factor=exp.(-τ*(norm(p)^2/(2m)-μ))
    normalized_data=zeros(diagram.max_order+1)

    for i in 1:diagram.max_order+1
        normalized_data[i]=data[i]*factor/zeroth_order
    end

    #data[:,index]=normalized_data

    return normalized_data
end

function jackknife(green_record,n_loop,ratio=0.1)


    num_bins=size(green_record[1],1)
    bin_width = Int(floor(ratio*n_loop)) #block binwidth
    jk_binwidth = n_loop-bin_width #jack knife binwidth, for complementary bins
    n_bins = cld(n_loop,bin_width)
    bin_variance=[]

    for i in 1:num_bins
        
        observables_array=[]
        for data in green_record
            append!(observables_array,data[i])
        end
        block_bins = collect(Iterators.partition(observables_array,bin_width)) 
        block_estimators = [sum(bin)/bin_width for bin in block_bins]
        jk_estimators = []

        observables_sum = sum(observables_array)
        for i in 1:n_bins
            jk_estimator = observables_sum - (bin_width * block_estimators[i])
            jk_estimator /= jk_binwidth
    
            append!(jk_estimators,jk_estimator)
        end
    
        #getting jack knife variance
        variance_jk = 0.0
        for k in 1:n_bins
            variance_jk = variance_jk+(jk_estimators[k] - mean(observables_array))^2
        end
        variance_jk *= (n_bins-1)/n_bins
        append!(bin_variance,variance_jk)
    end

    return bin_variance

end

function hist_measure!(diagram::Diagram,hist::Hist_Record,n_loop=5000,n_hist=10000,
                        p_ins=0.2,p_rem=0.2,p_swap=0.2,p_ex=0.2,p_to_0=0.2,p_from_0=0.5,p_ex_0=0.5)

    unnormalized_data=hist.unnormalized_data
    normalized_data=hist.normalized_data
    time_points=hist.time_points
    bin_width=hist.bin_width
    real_normalized=[p_ins,p_rem,p_swap,p_ex,p_to_0]
    real_normalized/=sum(real_normalized)
    fake_normalized=[p_from_0,p_ex_0]
    fake_normalized/=sum(fake_normalized)
    real_cumsum=cumsum(real_normalized)
    fake_cumsum=cumsum(fake_normalized)
    diagram.p_ins=real_normalized[1]
    diagram.p_rem=real_normalized[2]
    green_record=[]
    zero_record=[]
    regime=Diff_more()

    println("begin")
    for j in 1:n_loop
        println("loop.number:",j)
        for i in 1:n_hist
            q=rand()
            if diagram.order == 0
                if q<fake_cumsum[1]
                    diagram.p_ins=fake_normalized[1]
                    insert_arc!(diagram,regime)
                    diagram.p_ins=real_normalized[1]
                else
                    extend!(diagram)
                end
            else
                if q<real_cumsum[1]
                    insert_arc!(diagram,regime)
                elseif q<real_cumsum[2]
                    if diagram.order != 1
                        remove_arc!(diagram,regime) 
                    end       
                elseif q<real_cumsum[3]
                    swap_arc!(diagram)
                elseif q<real_cumsum[4]
                    extend!(diagram)
                else
                    if diagram.order == 1
                        diagram.p_ins=fake_normalized[1]
                        remove_arc!(diagram,regime)
                        diagram.p_ins=real_normalized[1]
                    end       
                end
            end
            unnormalized_data[diagram.order+1,Int(div(diagram.τ,bin_width,RoundUp))]+=1
        end

        
        # normalized_data=hist.normalized_data
        # green=normalized_data[1,:].*0
        # for i in 1:diagram.max_order+1
        #     green=green+normalized_data[i,:]
        # end
        green=unnormalized_data[1,:].*0
        for i in 1:diagram.max_order+1
            green=green+unnormalized_data[i,:]
        end
        push!(green_record,green)
        push!(zero_record,unnormalized_data[1,:])
    end

    hist.normalized_data=normalization(unnormalized_data,bin_width,diagram)
    normalized_data=hist.normalized_data
    bin_variance=jackknife(green_record,zero_record,n_loop,diagram,bin_width,0.1)

    return green_record,normalized_data,bin_variance
end

function hist_measure_2!(diagram::Diagram,hist::Hist_Record,n_loop=5000,n_hist=10000,
                        p_ins=0.2,p_rem=0.2,p_swap=0.2,p_to_0=0.2,p_from_0=1)

    unnormalized_data=hist.unnormalized_data
    normalized_data=hist.normalized_data
    time_points=hist.time_points
    bin_width=hist.bin_width
    real_normalized=[p_ins,p_rem,p_swap,p_to_0]
    real_normalized/=sum(real_normalized)
    fake_normalized=[p_from_0]
    fake_normalized/=sum(fake_normalized)
    real_cumsum=cumsum(real_normalized)
    fake_cumsum=cumsum(fake_normalized)
    diagram.p_ins=real_normalized[1]
    diagram.p_rem=real_normalized[2]
    green_record=[]
    zero_record=[]
    regime=Diff_more()

    println("begin")
    for j in 1:n_loop
        println("loop.number:",j)
        for i in 1:n_hist
            q=rand()
            if diagram.order == 0
                diagram.p_ins=fake_normalized[1]
                insert_arc!(diagram,regime)
                diagram.p_ins=real_normalized[1]
            else
                if q<real_cumsum[1]
                    insert_arc!(diagram,regime)
                elseif q<real_cumsum[2]
                    if diagram.order != 1
                        remove_arc!(diagram,regime) 
                    end       
                elseif q<real_cumsum[3]
                    swap_arc!(diagram)
                else
                    if diagram.order == 1
                        diagram.p_ins=fake_normalized[1]
                        remove_arc!(diagram,regime)
                        diagram.p_ins=real_normalized[1]
                    end       
                end
            end
            extend!(diagram)
            unnormalized_data[diagram.order+1,Int(div(diagram.τ,bin_width,RoundUp))]+=1
        end

        
        # normalized_data=hist.normalized_data
        # green=normalized_data[1,:].*0
        # for i in 1:diagram.max_order+1
        #     green=green+normalized_data[i,:]
        # end
        green=unnormalized_data[1,:].*0
        for i in 1:diagram.max_order+1
            green=green+unnormalized_data[i,:]
        end
        push!(green_record,green)
        push!(zero_record,unnormalized_data[1,:])
    end

    hist.normalized_data=normalization(unnormalized_data,bin_width,diagram)
    normalized_data=hist.normalized_data
    bin_variance=jackknife(green_record,zero_record,n_loop,diagram,bin_width,0.1)

    return green_record,normalized_data,bin_variance
end

function hist_measure_3!(diagram::Diagram,hist::Hist_Record,n_loop=5000,n_hist=10000,
                        p_ins=0.2,p_rem=0.2,p_to_0=0.2,p_from_0=1)

    unnormalized_data=hist.unnormalized_data
    normalized_data=hist.normalized_data
    time_points=hist.time_points
    bin_width=hist.bin_width
    real_normalized=[p_ins,p_rem,p_to_0]
    real_normalized/=sum(real_normalized)
    fake_normalized=[p_from_0]
    fake_normalized/=sum(fake_normalized)
    real_cumsum=cumsum(real_normalized)
    fake_cumsum=cumsum(fake_normalized)
    diagram.p_ins=real_normalized[1]
    diagram.p_rem=real_normalized[2]
    green_record=[]
    zero_record=[]
    regime=Diff_more()

    println("begin")
    for j in 1:n_loop
        println("loop.number:",j)
        for i in 1:n_hist
            q=rand()
            if diagram.order == 0
                diagram.p_ins=fake_normalized[1]
                insert_arc!(diagram,regime)
                diagram.p_ins=real_normalized[1]
            else
                if q<real_cumsum[1]
                    insert_arc!(diagram,regime)
                elseif q<real_cumsum[2]
                    if diagram.order != 1
                        remove_arc!(diagram,regime) 
                    end       
                else
                    if diagram.order == 1
                        diagram.p_ins=fake_normalized[1]
                        remove_arc!(diagram,regime)
                        diagram.p_ins=real_normalized[1]
                    end       
                end
            end
            extend!(diagram)
            unnormalized_data[diagram.order+1,Int(div(diagram.τ,bin_width,RoundUp))]+=1
        end

        
        # normalized_data=hist.normalized_data
        # green=normalized_data[1,:].*0
        # for i in 1:diagram.max_order+1
        #     green=green+normalized_data[i,:]
        # end
        green=unnormalized_data[1,:].*0
        for i in 1:diagram.max_order+1
            green=green+unnormalized_data[i,:]
        end
        push!(green_record,green)
        push!(zero_record,unnormalized_data[1,:])
    end

    hist.normalized_data=normalization(unnormalized_data,bin_width,diagram)
    normalized_data=hist.normalized_data
    bin_variance=jackknife(green_record,zero_record,n_loop,diagram,bin_width,0.1)

    return green_record,normalized_data,bin_variance
end

function discrete_measure!(diagram::Diagram,discrete::Discrete_Record,n_loop=1000,n_hist=2000,
                            p_ins=0.2,p_rem=0.2,p_swap=0.2,p_to_0=0.2,p_from_0=0.5)
    
    unnormalized_data=discrete.unnormalized_data
    normalized_data=discrete.normalized_data
    time_points=discrete.time_points
    real_normalized=[p_ins,p_rem,p_swap,p_to_0]
    real_normalized/=sum(real_normalized)
    fake_normalized=[p_from_0]
    fake_normalized/=sum(fake_normalized)
    real_cumsum=cumsum(real_normalized)
    fake_cumsum=cumsum(fake_normalized)
    diagram.p_ins=real_normalized[1]
    diagram.p_rem=real_normalized[2]
    green_record=[]
    regime=Diff_more()

    println("begin")
    for index in 1:length(time_points)

        diagram=Diagram(diagram.p[1], diagram.max_τ, diagram.max_order, diagram.mass, diagram.μ, diagram.ω, diagram.α)
        green_time_record=[]
        t=time_points[index]
        diagram.τ = t
        line_end=diagram.line_box[2*diagram.order+1]
        line_end.period[2]=t

        for j in 1:n_loop
            println("time_index:"*string(index)*" loop.number:"*string(j))
            for i in 1:n_hist
                q=rand()
                if diagram.order == 0
                    diagram.p_ins=fake_normalized[1]
                    insert_arc!(diagram,regime)
                    diagram.p_ins=real_normalized[1]
                else
                    if q<real_cumsum[1]
                        insert_arc!(diagram,regime)
                    elseif q<real_cumsum[2]
                        if diagram.order != 1
                            remove_arc!(diagram,regime) 
                        end       
                    elseif q<real_cumsum[3]
                        swap_arc!(diagram)
                    else
                        if diagram.order == 1
                            diagram.p_ins=fake_normalized[1]
                            remove_arc!(diagram,regime)
                            diagram.p_ins=real_normalized[1]
                        end
                    end
                end
                unnormalized_data[diagram.order+1,index]+=1
            end

            discrete.normalized_data[:,index]=normalization(unnormalized_data,time_points,index,diagram)
            normalized_data=discrete.normalized_data
            green=0

            for i in 1:diagram.max_order+1
                green=green+normalized_data[i,index]
            end

            push!(green_time_record,green)
        end
        push!(green_record,green_time_record)
    end


    # discrete.normalized_data=normalization(unnormalized_data,time_points,bin_width,diagram)
    # normalized_data=discrete.normalized_data
    # green=normalized_data[1,:].*0
    # for i in 1:diagram.max_order+1
    #     green=green+normalized_data[i,:]
    # end

    # push!(green_record,green)

    # bin_variance=jackknife(green_record,n_loop,0.1)

    return green_record,normalized_data#,bin_variance
end

function hist_measure2!(diagram::Diagram,hist::Hist_Record,folder, n_loop=5000, store_data=true, n_hist=100000,
                        p_ins=0.2,p_rem=0.2,p_from_0=1)

    if store_data
        address=joinpath(folder,"α="*string(round(diagram.α,digits=3)),
                        "k="*string(round(diagram.p[1],digits=3)),
                        "μ="*string(round(diagram.μ,digits=3)),
                        "histnum="*string(n_hist))
        mkpath(address)
        

        if isfile(joinpath(address,"diagram.jld2")) && isfile(joinpath(address,"hist.jld2"))
            diagram=load(joinpath(address,"diagram.jld2"), "diagram_a")
            hist=load(joinpath(address,"hist.jld2"), "hist_a")
            # @load joinpath(address,"diagram.jld2") diagram_a
            # @load joinpath(address,"hist.jld2") hist_a
        end
    end

    # diagram=diagram_a
    # hist=hist_a
    unnormalized_data=hist.unnormalized_data
    normalized_data=hist.normalized_data

    if store_data
        previous_data_0=copy(unnormalized_data[1,:])
        green=unnormalized_data[1,:].*0
        for i in 1:diagram.max_order+1
            green=green+unnormalized_data[i,:]
        end
        previous_data_t=green
    end

    time_points=hist.time_points
    bin_width=hist.bin_width
    real_normalized=[p_ins,p_rem]
    real_normalized/=sum(real_normalized)
    fake_normalized=[p_from_0]
    fake_normalized/=sum(fake_normalized)
    real_cumsum=cumsum(real_normalized)
    fake_cumsum=cumsum(fake_normalized)
    diagram.p_ins=real_normalized[1]
    diagram.p_rem=real_normalized[2]
    green_record=[]
    zero_record=[]
    regime=Diff_more()

    order=diagram.order
    m=diagram.mass
    μ=diagram.μ
    ω=diagram.ω
    α=diagram.α
    α_squared=2pi*α*sqrt(2)

    #println("begin")
    for j in 1:n_loop
        #println("loop.number:",j)
        for i in 1:n_hist
            q=rand()
            order=diagram.order
            if order == 0
                diagram.p_ins=fake_normalized[1]
                insert_arc!(diagram,regime)
                diagram.p_ins=real_normalized[1]
            else
                if q<real_cumsum[1]
                    insert_arc!(diagram,regime)      
                else
                    if order == 1
                        diagram.p_ins=fake_normalized[1]
                        remove_arc!(diagram,regime)
                        diagram.p_ins=real_normalized[1]
                    else
                        remove_arc!(diagram,regime)
                    end       
                end
            end
            extend!(diagram)
            unnormalized_data[order+1,Int(div(diagram.τ,bin_width,RoundUp))]+=1
        end

        green=unnormalized_data[1,:].*0
        for i in 1:diagram.max_order+1
            green=green+unnormalized_data[i,:]
        end

        if store_data
            CSV.write(joinpath(address,"total_green.csv"), DataFrame(transpose(hcat(green-previous_data_t)), :auto),append = true)
            previous_data_t=copy(green)
            CSV.write(joinpath(address,"zero_green.csv"), DataFrame(transpose(hcat(unnormalized_data[1,:]-previous_data_0)), :auto),append = true)
            previous_data_0=copy(unnormalized_data[1,:])
        end

        push!(green_record,green)
        push!(zero_record,unnormalized_data[1,:])
    end

    hist.normalized_data=normalization(unnormalized_data,bin_width,diagram)
    normalized_data=hist.normalized_data
    #bin_variance=jackknife(green_record,zero_record,n_loop,diagram,bin_width,0.1)

    if store_data
        save(joinpath(address,"diagram.jld2"), "diagram_a", diagram)
        save(joinpath(address,"hist.jld2"), "hist_a", hist)
    end
    # @save joinpath(address,"diagram.jld2") diagram_a=diagram
    # @save joinpath(address,"hist.jld2") hist_a=hist

    return diagram,hist,green_record,zero_record,normalized_data#,bin_variance#
end