include("Diagram.jl")
include("update.jl")
using Random
using CSV, DataFrames
using JLD2,FileIO

mutable struct Hist_Record

    bin_width:: Float64
    max_τ::Float64
    time_points::Array{Float64}
    max_order::Int64
    unnormalized_data::Array{Int64,2}
    normalized_data::Array{Float64,2}

    function Hist_Record(num_bins::Int64,max_τ::Real,max_order::Int64)
        bin_width=max_τ/num_bins
        unnormalized_data=zeros(max_order+1, num_bins)
        time_points=collect(1:num_bins)*bin_width.-(bin_width/2)
        new(bin_width, max_τ, time_points, max_order,unnormalized_data,unnormalized_data)
    end

end

function normalization(data,bin_width,diagram::Diagram)
    p=diagram.p
    μ=diagram.μ
    m=diagram.mass
    zeroth_order=deepcopy(data[1,:])
    # factor=exp.(-time_points.*(norm(p)^2/(2m)-μ))#/bin_width
    #factor=exp(-time_points[1]*(norm(p)^2/(2m)-μ))#/bin_width
    sum_zero=sum(zeroth_order)*bin_width
    factor=(1-exp(-diagram.max_τ*(norm(p)^2/(2m)-μ)))/(norm(p)^2/(2m)-μ)
    normalized_data=zeros(diagram.max_order+1,length(zeroth_order))

    for i in 1:diagram.max_order+1
        # normalized_data[i,:]=(data[i,:]./data[1,1]).*factor
        normalized_data[i,:]=(data[i,:]).*factor./sum_zero
    end

    return normalized_data
end

function binning(green_record,zero_record)
    
    total_green=copy(green_record)
    zero_green=copy(zero_record)

    for i in collect(length(total_green):-1:2)
        zero_green[i]-=zero_green[i-1]
        total_green[i]-=total_green[i-1]
    end

    return total_green, zero_green

end

function jackknife(green_record,zero_record,n_loop,diagram,binwidth,ratio=0.1)

    p=diagram.p
    μ=diagram.μ
    m=diagram.mass

    total_green, zero_green=binning(green_record,zero_record)
    num_bins=size(total_green[1],1)
    bin_width = Int(floor(ratio*n_loop)) #block binwidth
    jk_binwidth = n_loop-bin_width #jack knife binwidth, for complementary bins
    n_bins = cld(n_loop,bin_width)
    bin_variance=[]

    green_zero_array=[]
    for i in 1:n_loop
        append!(green_zero_array,sum(zero_green[i]))
    end
    block_bins_g0 = collect(Iterators.partition(green_zero_array,bin_width))
    block_estimators_g0 = [sum(bin) for bin in block_bins_g0]
    factor=(1-exp(-diagram.max_τ*(norm(p)^2/(2m)-μ)))/(norm(p)^2/(2m)-μ)
    sum_g0=sum(zero_record[n_loop])

    for i in 1:num_bins
        
        green_total_array=[]#total_green[:,i]
        # green_zero_array=[]#zero_green[:,i]
        for j in 1:n_loop
            append!(green_total_array,total_green[j][i])
            # append!(green_zero_array,zero_green[j][i])
        end

        block_bins_gt = collect(Iterators.partition(green_total_array,bin_width)) 
         
        block_estimators_gt = [sum(bin) for bin in block_bins_gt]
        
        jk_estimators = []

        sum_gt=green_record[n_loop][i]
        # sum_g0=zero_record[n_loop][i]

        # zeroth_order=deepcopy(data[1,:])
        # sum_zero=sum(zeroth_order)*bin_width
        
        # normalized_data=zeros(diagram.max_order+1,length(zeroth_order))
    
        # for i in 1:diagram.max_order+1
        #     # normalized_data[i,:]=(data[i,:]./data[1,1]).*factor
        #     normalized_data[i,:]=(data[i,:]).*factor./sum_zero
        # end
    
        # observables_sum = sum(observables_array)
        for j in 1:n_bins
            jk_estimator=(sum_gt-block_estimators_gt[j])/((sum_g0-block_estimators_g0[j])*binwidth)
            jk_estimator*=factor
            # jk_estimator = observables_sum - (bin_width * block_estimators[j])
            # jk_estimator /= jk_binwidth
            append!(jk_estimators,jk_estimator)
        end
        
        mean_value=sum_gt/(sum_g0*binwidth)*factor
        #getting jack knife variance
        variance_jk = 0.0
        for k in 1:n_bins
            variance_jk = variance_jk+(jk_estimators[k] - mean_value)^2
        end

        variance_jk *= (n_bins-1)/n_bins
        append!(bin_variance,variance_jk)
    end

    return bin_variance

end

function hist_measure!(diagram::Diagram,hist::Hist_Record,folder, n_loop=5000, store_data=true, n_hist=100000,
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
                else
                    if diagram.order == 1
                        diagram.p_ins=fake_normalized[1]
                        remove_arc!(diagram,regime)
                        diagram.p_ins=real_normalized[1]
                    else
                        remove_arc!(diagram,regime) 
                    end       
                end
            end
            extend!(diagram)
            unnormalized_data[diagram.order+1,Int(div(diagram.τ,bin_width,RoundUp))]+=1
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
    bin_variance=jackknife(green_record,zero_record,n_loop,diagram,bin_width,0.1)

    if store_data
        save(joinpath(address,"diagram.jld2"), "diagram_a", diagram)
        save(joinpath(address,"hist.jld2"), "hist_a", hist)
    end
    # @save joinpath(address,"diagram.jld2") diagram_a=diagram
    # @save joinpath(address,"hist.jld2") hist_a=hist

    return diagram,hist,green_record,zero_record,normalized_data,bin_variance#
end
