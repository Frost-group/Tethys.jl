include("Diagram.jl")
include("update.jl")
using Random
using LsqFit
using JLD2
# time update test
begin
    num_mea=1; regime=Diff_more(); regime_2=Diff_2()
    p_max=10; max_τ=30; max_order=500; mass=1; μ=-2.2; ω=1; α=2
    diagram_a=Diagram(p_max, max_τ, max_order, mass, μ, ω, α)

    loop=10
    record=[]
    for i in 1:loop
        τ_update!(diagram_a)
        append!(record,diagram_a.τ)
    end
    histogram(record)
end

# update functions test
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
        if remove_arc!(diagram_a,regime)
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


begin
    println("begin")
    # Random.seed!(1234)
    # record=[]
    loop=10000
    for j in 1:50000
        # loop_record=[]
        println(j)
        for i in 1:loop
            q=rand()
            if diagram_a.order == 0
                if q<0.5
                    diagram_a.p_rem=1/3
                    insert_arc!(diagram_a,regime)
                    diagram_a.p_rem=0.5
                else
                    extend!(diagram_a)
                end
            elseif diagram_a.order == 1
                if q<1/3
                    diagram_a.p_ins=1/3
                    diagram_a.p_rem=0.25
                    insert_arc!(diagram_a,regime)
                    diagram_a.p_ins=0.5
                    diagram_a.p_rem=0.5
                elseif q<2/3
                    diagram_a.p_rem=1/3
                    remove_arc!(diagram_a,regime)
                    diagram_a.p_rem=0.5
                else
                    extend!(diagram_a)
                end
            elseif  diagram_a.order == 2
                if q<0.25
                    insert_arc!(diagram_a,regime)
                elseif q<0.5
                    diagram_a.p_ins=1/3
                    diagram_a.p_rem=0.25
                    remove_arc!(diagram_a,regime)
                    diagram_a.p_ins=0.5
                    diagram_a.p_rem=0.5        
                elseif q<0.75
                    swap_arc!(diagram_a)
                else
                    extend!(diagram_a)
                end
            else
                if q<0.25
                    insert_arc!(diagram_a,regime)
                elseif q<0.5
                    remove_arc!(diagram_a,regime)        
                elseif q<0.75
                    swap_arc!(diagram_a)
                else
                    extend!(diagram_a)
                end
            end
            statis[diagram_a.order+1,Int(div(diagram_a.τ,bin_width,RoundUp))]+=1
        end
    end
end

# measurement test (with separating order=0,1 cases)
begin
    println("begin")
    # Random.seed!(1234)
    # record=[]
    diagram_a=Diagram(0, max_τ, max_order, mass, μ, ω, α)
    diagram_a.p=[0.,0.,0.]
    loop=10000
    for j in 1:5000
        # loop_record=[]
        println(j)
        for i in 1:loop
            q=rand()
            if diagram_a.order == 0
                # println("0 insert")
                insert_arc!(diagram_a,regime)
                extend!(diagram_a)
                # if q<0.5
                #     # println("0insert")
                #     insert_arc!(diagram_a,regime)
                # else
                #     extend!(diagram_a)
                # end
            elseif diagram_a.order == 1
                if q<0.5
                    # println("insert")
                    insert_arc!(diagram_a,regime)
                    extend!(diagram_a)
                # elseif q<2/3
                #     # println("remove")
                #     remove_arc!(diagram_a,regime)
                else
                    # println("remove")
                    remove_arc!(diagram_a,regime)
                    extend!(diagram_a)
                end
            else
                if q<0.5
                    insert_arc!(diagram_a,regime)
                    swap_arc!(diagram_a)
                    extend!(diagram_a)
                else#if q<0.5
                    remove_arc!(diagram_a,regime)
                    swap_arc!(diagram_a)
                    extend!(diagram_a)
                # elseif q<0.75
                #     swap_arc!(diagram_a)
                # else
                #     # swap_arc!(diagram_a)
                #     extend!(diagram_a)
                end

            end
            if j>10
                statis[diagram_a.order+1,Int(div(diagram_a.τ,bin_width,RoundUp))]+=1
            end
            # push!(loop_record,measure(diagram_a))

            # if i >=100
            #     push!(loop_record,measure(diagram_a))
            # end
            # if diagram_a.order==0
            #     println("zero")
            #     # println(diagram_a.order)
            # end
        end
        # append!(record,loop_record)
    end
    # save_object("E://data_record//record"*string(num_mea)*".jld2", record)
    # num_mea+=1
end

# measurement test (without separating order=0,1 cases)
begin
    regime_2=Diff_2()
    println("begin")
    record=[]
    loop=10000
    for j in 1:1000
        diagram_a=Diagram(p_max, max_τ, max_order, mass, μ, ω, α)
        diagram_a.p=[0.,0.,0.]
        loop_record=[]
        println(j)
        for i in 1:loop
            q=rand()
            if q<0.25
                insert_arc!(diagram_a)
            elseif q<0.5
                remove_arc!(diagram_a,regime_2)
            elseif q<0.75
                swap_arc!(diagram_a)
            else
                extend!(diagram_a)
            end
            push!(loop_record,measure(diagram_a))
        end
        append!(record,loop_record)
    end
end

# tau histogram for order=0-10
begin
    bins=range(0, stop = 5, length = 300)
    mod_record=[box[2] for box in record if box[1]==0]
    histogram(mod_record,bins=bins,label="order="*string(0))
    max_order=10
    for i in 1:max_order
        mod_record=[box[2] for box in record if box[1]==i]
        if i == max_order
            display(histogram!(mod_record,bins=bins,label="order="*string(i),xlabel="τ",ylabel="number of data"))
        else
            histogram!(mod_record,bins=bins,label="order="*string(i))
        end
    end
    # savefig("different_order.png")
end

# tau histogram for order=0
begin
    bins=range(0, stop = 5, length = 100)
    mod_record=[box[2] for box in record if box[1]==0]
    histogram(mod_record,bins=bins,label="order="*string(0),xlabel="τ",ylabel="number of data")
    # savefig("zero_order.png")
end

# tau histogram for all the orders
begin
    bins=range(0, stop = 40, length = 300)
    mod_record=[box[2] for box in record]
    histogram(mod_record,bins=bins,yaxis=:log,label="all the orders",xlabel="τ",ylabel="number of data")
    # savefig("all_order.png")
end

# order histogram
begin
    mod_record=[box[1] for box in record]
    histogram(mod_record,yaxis=:log,xlabel="order",ylabel="number of data")
    # savefig("order_hist.png")
end

begin
    record=[2,3,4,4,3]
    save_object("record1.jld2", record)
    b=load_object("E://data_record//record"*string(num_mea)*".jld2")
    # save_object()
end