include("Diagram.jl")

function insert_arc!(diagram::Diagram,regime::Diff_2)

    order=diagram.order
    m=diagram.mass
    μ=diagram.μ
    ω=diagram.ω
    α=diagram.α
    α_squared=2pi*α*sqrt(2)

    if order+1>diagram.max_order
        return false
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
    new_arc=Arc(q,[τ_1,τ_2],ω,index,index+2)

    time=[τ_L,τ_1,τ_2,τ_R]
    k_box=[line.k,line.k-q,line.k]
    covered=[false,true,false]
    line_to_add=[]

    w_x=green_zero(line)
    w_y=phonon_propagator(new_arc)*α_squared/(2*pi)^3
    # println("w_x",w_x)
    # println("w_y",w_y)
    for i in 1:3
        new_line=Line(k_box[i] ,[time[i],time[i+1]], m, μ, index+i-1, covered[i])
        push!(line_to_add,new_line)
        w_y*=green_zero(new_line)
    end

    p_x_y=diagram.p_ins/(τ_R-τ_L)*ω*exp(-ω*arc_T)
    p_x_y*=exp(-norm(q)^2/(2m)*arc_T)/(2pi*m/arc_T)^1.5

    p_y_x=diagram.p_rem/(order+1)
    
    r=w_y*p_y_x/(w_x*p_x_y)*(2order+1)#*(1-exp(-ω*(τ_R-τ_1)))
    
    # vertex=4*pi*α*ω^1.5/sqrt(2m)
    # r=sqrt(m/(2*pi*arc_T))^3*vertex*(2order+1)*(τ_R-τ_L)*diagram.p_rem/diagram.p_ins*exp(arc_T*(dot(q,line.k)/m+ω))/(ω*(order+1)*norm(q)^2)
    # println("insert_r=",r)

    if r<rand()
        return false
    else
        diagram.order+=1
        deleteat!(line_box, index)
        sign_box=diagram.sign_box
        sign_to_add=[[copy(sign_box[index][1]),-1],[-1,1],[1,copy(sign_box[index][2])]]
        deleteat!(sign_box, index)

        for i in 1:3
            insert!(line_box, index, line_to_add[4-i])
            insert!(sign_box, index, sign_to_add[4-i])
        end
        # println("insert index is ",[index,index+2])
        if length(line_box)>=index+3
            for i in index+3:length(line_box)
                line_box[i].index=i
            end
        end

        for arc in diagram.arc_box
            if arc.index_out<=index
                continue
            elseif arc.index_in>=index
                arc.index_in+=2
                arc.index_out+=2
            elseif arc.index_in<index && arc.index_out>index
                arc.index_out+=2
            end
        end

        push!(diagram.arc_box,new_arc)

        return true
    end
end

function remove_arc!(diagram::Diagram,regime::Diff_2)

    order=diagram.order
    m=diagram.mass
    μ=diagram.μ
    ω=diagram.ω
    α=diagram.α
    α_squared=2pi*α*sqrt(2)

    if order-1<0
        return false
    end
    
    arc_box=diagram.arc_box
    index=rand(1:length(arc_box))
    arc=arc_box[index]
    index_in=arc.index_in
    index_out=arc.index_out
    q=arc.q

    if index_out-index_in>2
        return false
    end

    line_box=diagram.line_box
    line_in=line_box[index_in]
    line_out=line_box[index_out]
    line_to_rem=[line_box[index_in],line_box[index_in+1],line_box[index_in+2]]

    τ_L=line_in.period[1]
    τ_R=line_out.period[2]
    new_line=Line(line_in.k ,[τ_L,τ_R], m, μ, index_in,false)

    w_x=green_zero(new_line)
    w_y=phonon_propagator(arc)*α_squared/(2*pi)^3

    for i in 1:3
        w_y*=green_zero(line_to_rem[i])
    end

    τ_1=copy(arc.period[1])
    τ_2=copy(arc.period[2])
    arc_T=τ_2-τ_1
    p_x_y=diagram.p_ins/(τ_R-τ_L)*ω*exp(-ω*arc_T)
    p_x_y*=exp(-norm(q)^2/(2m)*arc_T)/(2pi*m/arc_T)^1.5

    p_y_x=diagram.p_rem/order

    r=(w_x*p_x_y)/(w_y*p_y_x)/(2order-1)
    # println("remove_r=",r)
    # vertex=4*pi*α*ω^1.5/sqrt(2m)
    # r=sqrt(m/(2*pi*arc_T))^3*vertex*(2order-1)*(τ_R-τ_L)*diagram.p_rem/diagram.p_ins*exp(arc_T*(dot(q,line_in.k)/m+ω))/(ω*(order)*norm(q)^2)
    # r=1/r

    if r<rand()
        return false
    else
        diagram.order-=1
        
        deleteat!(arc_box, index)
        deleteat!(line_box, index_in:index_out)
        insert!(line_box, index_in, new_line)
        covered=false

        sign_box=diagram.sign_box
        sign_to_add=[sign_box[index_in][1],sign_box[index_out][2]]
        deleteat!(sign_box, index_in:index_out)
        insert!(sign_box, index_in, sign_to_add)
        

        if length(line_box)>=index_in+1
            for i in index_in+1:length(line_box)
                line_box[i].index=i
            end
        end

        for arc in diagram.arc_box
            if arc.index_out<=index_in
                continue
            elseif arc.index_in>=(index_in+2)
                arc.index_in-=2
                arc.index_out-=2
            elseif arc.index_in<index_in && arc.index_out>index_in
                arc.index_out-=2
                covered=true
            # elseif arc.index_in>index_in && arc.index_out<index_out
            #     arc.index_in-=1
            #     arc.index_out-=1
            end
        end

        if covered
            line_box[index_in].covered=covered
        end

        return true
    end
end



