include("Diagram.jl")
include("update.jl")
using Random

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

function normalization(data,time_points,bin_width,diagram::Diagram)
    p=diagram.p
    μ=diagram.μ
    m=diagram.mass
    zeroth_order=deepcopy(data[1,:])
    # factor=exp.(-time_points.*(norm(p)^2/(2m)-μ))#/bin_width
    factor=exp(-time_points[1]*(norm(p)^2/(2m)-μ))#/bin_width
    normalized_data=zeros(diagram.max_order+1,length(zeroth_order))

    for i in 1:diagram.max_order+1
        normalized_data[i,:]=(data[i,:]./data[1,1]).*factor
    end

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

        hist.normalized_data=normalization(unnormalized_data,time_points,bin_width,diagram)
        normalized_data=hist.normalized_data
        green=normalized_data[1,:].*0
        for i in 1:diagram.max_order+1
            green=green+normalized_data[i,:]
        end
        push!(green_record,green)
    end

    bin_variance=jackknife(green_record,n_loop,0.1)

    return green_record,normalized_data,bin_variance
end