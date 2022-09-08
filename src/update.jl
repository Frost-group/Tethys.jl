include("Diagram.jl")

function p_update!(diagram::Diagram)

    if diagram.order != 0
        return 
    end

    p_sample=diagram.p_grid
    p_new=[rand(p_sample)*rand([-1,1]),0,0]
    p_old=copy(diagram.p)
    old_green=green_zero(diagram)
    diagram.p=p_new
    new_green=green_zero(diagram)

    if new_green/old_green >= rand()
        return true
    else
        diagram.p=p_old
        return false
    end
end

function τ_update!(diagram::Diagram)

    if diagram.order != 0
        return false
    end

    p=diagram.p
    μ=diagram.μ
    m=diagram.mass
    dispersion=norm(p)^2/(2m)-μ

    τ_new=-log(rand())/abs(dispersion)

    if τ_new>=diagram.max_τ
        return false
    else
        diagram.τ=τ_new
        diagram.line_box[1].period[2]=τ_new
        return true
    end
end

function insert_arc!(diagram::Diagram)

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

function insert_arc!(diagram::Diagram,regime::Diff_more)

    order=diagram.order
    m=diagram.mass
    μ=diagram.μ
    ω=diagram.ω
    α=diagram.α
    α_squared=2pi*α*sqrt(2)
    τ=diagram.τ

    if order+1>diagram.max_order
        return false
    end

    line_box=diagram.line_box
    index_in=rand(1:length(line_box))
    line=line_box[index_in]
    τ_L=deepcopy(line.period[1])
    τ_R=deepcopy(line.period[2])
    k_in=deepcopy(line.k)

    τ_1=rand(Uniform(τ_L,τ_R))
    τ_2=τ_1-log(rand())/ω 

    if τ_2 > τ
        return false
    end

    line.period[1]=τ_1
    arc_T=τ_2-τ_1
    q=rand(Normal(0,sqrt(m/arc_T)),3)
    
    w_x=1
    w_y=1
    τ_R_2=0
    index_out=0
    k_out=0

    #not set covered yet
    for i in index_in:2order+1
        line_tem=line_box[i]
        if line_tem.period[2]<τ_2
            w_x*=green_zero(line_tem)
            line_tem.k-=q
            line_tem.index+=1
            w_x/=green_zero(line_tem)
            continue
        else
            k_out=deepcopy(line_tem.k)
            τ_R_2=deepcopy(line_tem.period[2])
            line_tem.period[2]=τ_2
            w_x*=green_zero(line_tem)
            line_tem.k-=q
            line_tem.index+=1
            w_x/=green_zero(line_tem)
            index_out=i+2
            break
        end
    end

    # println("index_in:",index_in)
    # println("index_out",index_out-2)
    new_arc=Arc(q,[τ_1,τ_2],ω,index_in,index_out)

    w_y*=phonon_propagator(new_arc)*α_squared/(2*pi)^3
    p_x_y=diagram.p_ins/(2order+1)/(τ_R-τ_L)*ω*exp(-ω*arc_T)#/(1-exp(-ω*(τ-τ_1)))#
    p_x_y*=exp(-norm(q)^2/(2m)*arc_T)/(2pi*m/arc_T)^1.5#/normalization(τ,τ_R,τ_L,ω)#*weighting/exp(-ω*τ_1)
    # println("ratio=",line_length/(τ-τ_1))
    p_y_x=diagram.p_rem/(order+1)
    r=w_y*p_y_x/(w_x*p_x_y)#^2
    # println("normal",normalization(τ,τ_R,τ_L,ω))
    # println("insert_r=",r)
    

    if r<rand()
        line_box[index_in].period[1]=τ_L
        line_box[index_out-2].period[2]=τ_R_2
        for i in index_in:index_out-2
            line_tem=line_box[i]
            line_tem.index=i
            line_tem.k+=q
        end

        return false
    else
        diagram.order+=1

        if index_out-index_in==2
            line_box[index_in].covered=true
        else
            line_box[index_in].covered=false
            line_box[index_out-2].covered=false
        end

        line_tem=Line(k_in,[τ_L,τ_1], m, μ, index_in, false)
        insert!(line_box, index_in, line_tem)
        line_tem=Line(k_out,[τ_2,τ_R_2], m, μ, index_out, false)
        insert!(line_box, index_out, line_tem)

        sign_box=diagram.sign_box
        if index_out-index_in==2
            sign_to_add=[[copy(sign_box[index_in][1]),-1],[-1,1],[1,copy(sign_box[index_in][2])]]
            deleteat!(sign_box, index_in)
            for i in 1:3
                insert!(sign_box, index_in, sign_to_add[4-i])
            end
            
        else

            sign_to_add=[[copy(sign_box[index_out-2][1]),1],[1,copy(sign_box[index_out-2][2])]]
            deleteat!(sign_box, index_out-2)
            for i in 1:2
                insert!(sign_box, index_out-2, sign_to_add[3-i])
            end

            sign_to_add=[[copy(sign_box[index_in][1]),-1],[-1,copy(sign_box[index_in][2])]]
            deleteat!(sign_box, index_in)
            for i in 1:2
                insert!(sign_box, index_in, sign_to_add[3-i])
            end
        end

        if length(line_box)>=index_out+1
            for i in index_out+1:length(line_box)
                line_box[i].index=i
            end
        end

        for arc in diagram.arc_box
            if arc.index_out<=index_in
                continue
            elseif arc.index_in>=index_out-2
                arc.index_in+=2
                arc.index_out+=2
            elseif arc.index_in<index_in && arc.index_out>index_out-2
                arc.index_out+=2
            elseif arc.index_in>=index_in && arc.index_out<=index_out-2
                arc.index_in+=1
                arc.index_out+=1
            elseif index_in<=arc.index_in<index_out-2
                arc.index_in+=1
                arc.index_out+=2
            elseif index_in<arc.index_out<=index_out-2
                arc.index_out+=1
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

function remove_arc!(diagram::Diagram,regime::Diff_more)

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

    # if index_out-index_in>2
    #     return false
    # end

    # old_box=deepcopy(diagram.line_box)
    line_box=diagram.line_box
    line_in=line_box[index_in]
    line_out=line_box[index_out]
    # line_to_rem=[line_box[index_in],line_box[index_in+1],line_box[index_in+2]]

    τ_L=line_in.period[1]

    if index_out-index_in==2
        τ_R=line_out.period[2]
    else
        τ_R=line_box[index_in+1].period[2]
    end

    τ_R_2=line_out.period[2]
    # new_line=Line(line_in.k ,[τ_L,τ_R], m, μ, index_in,false)
    
    w_x=1
    w_y=phonon_propagator(arc)*α_squared/(2*pi)^3

    for i in index_in+1:index_out-1
        line_tem=line_box[i]
        w_y*=green_zero(line_tem)
        line_tem.index-=1
        line_tem.k+=q
        w_y/=green_zero(line_tem)
    end

    τ_1=arc.period[1]
    τ_2=arc.period[2]
    arc_T=τ_2-τ_1
    p_x_y=diagram.p_ins/(2order-1)/(τ_R-τ_L)*ω*exp(-ω*arc_T)#/(1-exp(-ω*(diagram.τ-τ_1)))#/(2order-1)*(τ_R_2-line_box[index_out-1].period[1])/(diagram.τ-τ_1)#diagram.τ#
    p_x_y*=exp(-norm(q)^2/(2m)*arc_T)/(2pi*m/arc_T)^1.5#/normalization(diagram.τ,τ_R,τ_L,ω)
    # p_x_y*=exp(-ω*line_box[index_out-1].period[1])-exp(-ω*τ_R_2)
    # p_x_y/=exp(-ω*τ_1)

    p_y_x=diagram.p_rem/order

    r=(w_x*p_x_y)/(w_y*p_y_x)#^2
    # println("remove_r=",r)
    # vertex=4*pi*α*ω^1.5/sqrt(2m)
    # r=sqrt(m/(2*pi*arc_T))^3*vertex*(2order-1)*(τ_R-τ_L)*diagram.p_rem/diagram.p_ins*exp(arc_T*dot(q,line_in.k)/m)/(ω*(order)*norm(q)^2)
    # r=1/r

    if r<rand()
        for i in index_in+1:index_out-1
            line_tem=line_box[i]
            line_tem.index+=1
            line_tem.k-=q
        end
        # diagram.line_box=old_box
        return false
    else
        sign_box=diagram.sign_box
        diagram.order-=1
        deleteat!(arc_box, index)

        if index_out-index_in==2
            line_tem=line_box[index_in+1]
            line_tem.period[1]=τ_L
            line_tem.period[2]=τ_R_2
            line_tem.covered=false
            sign_to_add=[sign_box[index_in][1],sign_box[index_out][2]]
            deleteat!(sign_box, index_in:index_out)
            insert!(sign_box, index_in, sign_to_add)
            # line_box[index_in+1]=line_tem
        else
            line_tem=line_box[index_in+1]
            line_tem.period[1]=τ_L
            sign_to_add=[sign_box[index_in][1],sign_box[index_in+1][2]]
            deleteat!(sign_box, index_in:index_in+1)
            insert!(sign_box, index_in, sign_to_add)
            # line_box[index_in+1]=line_tem

            line_tem=line_box[index_out-1]
            line_tem.period[2]=τ_R_2
            sign_to_add=[sign_box[index_out-2][1],sign_box[index_out-1][2]]
            deleteat!(sign_box, index_out-2:index_out-1)
            insert!(sign_box, index_out-2, sign_to_add)
            # line_box[index_in+1]=line_tem
        end

        deleteat!(line_box, [index_in,index_out])

        # deleteat!(line_box, index_in:index_out)
        # insert!(line_box, index_in, new_line)
        # covered=false

        # sign_box=diagram.sign_box
        # sign_to_add=[sign_box[index_in][1],sign_box[index_out][2]]
        # deleteat!(sign_box, index_in:index_out)
        # insert!(sign_box, index_in, sign_to_add)
        

        if length(line_box)>=index_out+1
            for i in index_out+1:length(line_box)
                line_box[i].index=i
            end
        end

        for arc in diagram.arc_box
            if arc.index_out<=index_in
                continue
            elseif arc.index_in>=index_out
                arc.index_in-=2
                arc.index_out-=2
            elseif arc.index_in<index_in && arc.index_out>index_out
                arc.index_out-=2
                if arc.index_out-arc.index_in == 2
                    line_box[arc.index_in+1].covered=true
                end
            elseif arc.index_in>index_in && arc.index_out<index_out
                arc.index_in-=1
                arc.index_out-=1
            elseif index_in<arc.index_in<index_out
                arc.index_in-=1
                arc.index_out-=2
                if arc.index_out-arc.index_in == 2
                    line_box[arc.index_in+1].covered=true
                end
            elseif index_in<arc.index_out<index_out
                arc.index_out-=1
                if arc.index_out-arc.index_in == 2
                    line_box[arc.index_in+1].covered=true
                end
                # covered=true
            # elseif arc.index_in>index_in && arc.index_out<index_out
            #     arc.index_in-=1
            #     arc.index_out-=1
            end
        end

        return true
    end
end

function arc_judge(arc::Arc,sign::Int64,bound::Bool,index::Int64)
    # bound true is right, false is left
    if bound 
        if sign == 1
            return arc.index_out-1 == index
        else
            return arc.index_in == index
        end
    else
        if sign == 1
            return arc.index_out == index
        else
            return arc.index_in+1 == index
        end
    end
end

function swap_arc!(diagram::Diagram)

    order=diagram.order
    m=diagram.mass
    μ=diagram.μ
    ω=diagram.ω

    if order<2
        return false
    end

    line_box=diagram.line_box
    line_index=rand(2:2*order)
    chosen_line=line_box[line_index]

    if chosen_line.covered
        return false
    end

    sign=diagram.sign_box[line_index]
    arc_box=diagram.arc_box
    left_check=false
    right_check=false
    left_index=0
    right_index=0

    for i in 1:order
        arc=arc_box[i]
        if !left_check
            if arc_judge(arc,sign[1],false,line_index)
                left_index=i
                left_check=true
            end
        end

        if !right_check
            if arc_judge(arc,sign[2],true,line_index)
                right_index=i
                right_check=true
            end
        end

        if right_check && left_check
            break
        end
    end
    # println("swap_index is:",line_index)
    arc_l=arc_box[left_index]
    arc_r=arc_box[right_index]

    q1=arc_l.q
    q2=arc_r.q

    if sign[1] == 1
        new_arc_l=Arc(q1,[arc_l.period[1],chosen_line.period[2]],ω,arc_l.index_in,line_index+1)
    else
        new_arc_l=Arc(q1,[chosen_line.period[2],arc_l.period[2]],ω,line_index,arc_l.index_out)
    end

    if sign[2] == 1
        new_arc_r=Arc(q2,[arc_r.period[1],chosen_line.period[1]],ω,arc_r.index_in,line_index)
    else
        new_arc_r=Arc(q2,[chosen_line.period[1],arc_r.period[2]],ω,line_index-1,arc_r.index_out)
    end

    new_line=Line(chosen_line.k-sign[1]*q1+sign[2]*q2 ,chosen_line.period, m, μ, line_index,false)
    w_x=green_zero(chosen_line)*phonon_propagator(arc_l)*phonon_propagator(arc_r)
    w_y=green_zero(new_line)*phonon_propagator(new_arc_l)*phonon_propagator(new_arc_r)

    r=w_y/w_x

    if r<rand()
        return false
    else
        deleteat!(line_box, line_index)
        insert!(line_box, line_index, new_line)
        sign_box=diagram.sign_box
        deleteat!(sign_box, line_index)
        insert!(sign_box, line_index, [sign[2],sign[1]])
        sign_box[line_index-1]=[sign_box[line_index-1][1],sign[2]]
        sign_box[line_index+1]=[sign[1],sign_box[line_index+1][2]]

        deleteat!(arc_box, left_index)
        insert!(arc_box, left_index, new_arc_l)
        deleteat!(arc_box, right_index)
        insert!(arc_box, right_index, new_arc_r)

        if new_arc_l.index_out-new_arc_l.index_in == 2
            line_box[new_arc_l.index_in+1].covered=true
        end

        if new_arc_r.index_out-new_arc_r.index_in == 2
            line_box[new_arc_r.index_in+1].covered=true
        end
        return true
    end
end

function extend!(diagram::Diagram)

    p=diagram.p
    μ=diagram.μ
    m=diagram.mass
    dispersion=norm(p)^2/(2m)-μ
    line_box=diagram.line_box
    order=diagram.order
    line_end=line_box[2*order+1]

    τ_new=line_end.period[1]-log(rand())/abs(dispersion)

    if τ_new>diagram.max_τ
        return false
    end

    line_end.period[2]=τ_new
    line_box[2*order+1]=line_end
    diagram.τ=τ_new

    return true

end
