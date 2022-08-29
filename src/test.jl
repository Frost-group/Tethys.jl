include("Diagram.jl")

# time update test
begin
    p_max=10; max_τ=40; max_order=500; mass=1; μ=-6; ω=1; α=5
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

# measurement test (with separating order=0,1 cases)
begin
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
                    remove_arc!(diagram_a)
                else
                    extend!(diagram_a)
                end
            else
                if q<0.25
                    insert_arc!(diagram_a)
                elseif q<0.5
                    remove_arc!(diagram_a)
                elseif q<0.75
                    swap_arc!(diagram_a)
                else
                    extend!(diagram_a)
                end
            end
            push!(loop_record,measure(diagram_a))

            # if i >=100
            #     push!(loop_record,measure(diagram_a))
            # end
            # if diagram_a.order==0
            #     println("zero")
            #     # println(diagram_a.order)
            # end
        end
        append!(record,loop_record)
    end
end

# measurement test (without separating order=0,1 cases)
begin
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
                remove_arc!(diagram_a)
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
    bins=range(0, stop = 5, length = 100)
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
end

# tau histogram for order=0
begin
    bins=range(0, stop = 5, length = 100)
    mod_record=[box[2] for box in record if box[1]==0]
    histogram(mod_record,bins=bins,label="order="*string(0),xlabel="τ",ylabel="number of data")
end

# tau histogram for all the orders
begin
    bins=range(0, stop = 40, length = 50)
    mod_record=[box[2] for box in record]
    histogram(mod_record,bins=bins,yaxis=:log,xlabel="τ",ylabel="number of data")

end

# order histogram
begin
    mod_record=[box[1] for box in record]
    histogram(mod_record,yaxis=:log,xlabel="order",ylabel="number of data")
end
