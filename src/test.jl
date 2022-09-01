include("Diagram.jl")
using Random
using JLD2
# time update test
begin
    num_mea=1; regime=Diff_more()
    p_max=10; max_τ=40; max_order=0; mass=1; μ=-6; ω=1; α=5
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
    statis=zeros(max_order+1, 300)
    bin_width=max_τ/300
end

# measurement test (with separating order=0,1 cases)
begin
    println("begin")
    # Random.seed!(124)
    # record=[]
    loop=1000
    for j in 1:100000
        diagram_a=Diagram(p_max, max_τ, max_order, mass, μ, ω, α)
        diagram_a.p=[0.,0.,0.]
        # loop_record=[]
        println(j)
        for i in 1:loop
            # if j==9
            #     println(diagram_a.order)
            # end
            # println("j:",j,",i:",i)
            q=rand()
            if diagram_a.order == 0
                if q<0.5
                    insert_arc!(diagram_a)
                else
                    extend!(diagram_a)
                end
            elseif diagram_a.order == 1
                if q<1/3
                    insert_arc!(diagram_a)
                elseif q<2/3
                    remove_arc!(diagram_a,regime)
                else
                    extend!(diagram_a)
                end
            else
                if q<0.25
                    insert_arc!(diagram_a)
                elseif q<0.5
                    remove_arc!(diagram_a,regime)
                elseif q<0.75
                    swap_arc!(diagram_a)
                else
                    extend!(diagram_a)
                end

            end
            statis[diagram_a.order+1,Int(div(diagram_a.τ,bin_width,RoundUp))]+=1
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

begin
    time=collect(1:75)*bin_width.-(bin_width/2)
    plot(time,statis[1,1:75],yaxis=:log)
    display(plot!(time,statis[1,1]*exp.(μ.*time)))
    println(-(log(statis[1,1])- log(statis[1,20]))/(bin_width*19))
    println(-(log(statis[1,75])- log(statis[1,150]))/(bin_width*74))
end

begin
    time=collect(1:100)*bin_width.-(bin_width/2)
    plot(time,statis[1,1:100])
    for order in 2:10+1
        if order == 10+1
            display(plot!(time,statis[order,1:100]))
        else
            plot!(time,statis[order,1:100])
        end
    end
end

begin
    max_check=20
    data=[sum(statis[i,:]) for i in 1:max_check+1]
    plot(collect(0:max_check),data, markershape=:circle,yaxis=:log,xticks=0:max_check)
end

begin
    time=collect(1:300)*bin_width.-(bin_width/2)
    data=[sum(statis[:,i]) for i in 1:300]
    display(plot(time,data,markershape=:circle,yaxis=:log))
    println(-(log(data[150])- log(data[75]))/(bin_width*74)+μ)
end