include("Diagram.jl")

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


function insert_arc!(diagram::Diagram,order::Int64,m::Float64,μ::Float64,ω::Float64,α_squared::Float64)

    τ=diagram.τ
    cross_over=false

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

    arc_T=τ_2-τ_1

    if arc_T > τ
        return false
    end

    line.period[1]=τ_1
    arc_T=τ_2-τ_1
    q=MVector{3}(rand(Normal(0,sqrt(m/arc_T)),3))
    
    w_x=1.0
    w_y=1.0
    total_dis=0.0
    τ_R_2=0.0
    index_out=0
    k_out=0

    if τ_2<τ
        #not set covered yet
        for i in index_in:2order+1
            line_tem=line_box[i]
            # total_dis+=dispersion(line_tem)
            # line_tem.k-=q
            # line_tem.index+=1
            # total_dis-=dispersion(line_tem)
            if line_tem.period[2]<τ_2
                total_dis+=dispersion(line_tem)
                line_tem.k-=q
                line_tem.index+=1
                total_dis-=dispersion(line_tem)
                continue
            else
                k_out=deepcopy(line_tem.k)
                τ_R_2=deepcopy(line_tem.period[2])
                line_tem.period[2]=τ_2
                total_dis+=dispersion(line_tem)
                line_tem.k-=q
                line_tem.index+=1
                total_dis-=dispersion(line_tem)
                index_out=i+2
                break
            end
        end
    else
        cross_over=true
        τ_2=τ_2-τ
        for i in index_in:2order+1
            line_tem=line_box[i]
            total_dis+=dispersion(line_tem)
            line_tem.k-=q
            line_tem.index+=1
            total_dis-=dispersion(line_tem)
        end

        for i in 1:index_in
            line_tem=line_box[i]
            # total_dis+=dispersion(line_tem)
            # line_tem.k-=q
            # total_dis-=dispersion(line_tem)

            if line_tem.period[2]<τ_2
                total_dis+=dispersion(line_tem)
                line_tem.k-=q
                total_dis-=dispersion(line_tem)
                continue
            elseif i != index_in
                k_out=deepcopy(line_tem.k)
                τ_R_2=deepcopy(line_tem.period[2])
                line_tem.period[2]=τ_2
                total_dis+=dispersion(line_tem)
                line_tem.k-=q
                total_dis-=dispersion(line_tem)
                index_out=i+1
                break
            else
                k_out=deepcopy(line_tem.k)
                line_tem.period[1]=τ_L
                line_tem.period[2]=τ_2
                total_dis-=dispersion(line_tem)
                line_tem.k+=q
                total_dis+=dispersion(line_tem)
                index_out=i+1
                line_tem.period[1]=τ_2
                line_tem.period[2]=τ_1
            end
        end

        index_in=index_in+1
    end


    new_arc=Arc(q,[τ_1,τ_2],ω,index_in,index_out)

    p_x_y=diagram.p_ins/(2order+1)/(τ_R-τ_L)
    p_x_y*=exp(-norm(q)^2/(2m)*arc_T)/(2pi*m/arc_T)^1.5
    p_y_x=diagram.p_rem/(order+1)
    r=α_squared*p_y_x/(exp(total_dis)*p_x_y*(2*pi)^3*norm(q)^2)#*(1-exp(-ω*τ))/ω
    
    if r<rand()
        if !cross_over
            line_box[index_in].period[1]=τ_L
            line_box[index_out-2].period[2]=τ_R_2
            for i in index_in:index_out-2
                line_tem=line_box[i]
                line_tem.index=i
                line_tem.k+=q
            end
        else
            if index_in!=index_out
                line_box[index_in-1].period[1]=τ_L
                line_box[index_out-1].period[2]=τ_R_2
                for i in index_in-1:2order+1
                    line_tem=line_box[i]
                    line_tem.index=i
                    line_tem.k+=q
                end

                for i in 1:index_out-1
                    line_tem=line_box[i]
                    line_tem.index=i
                    line_tem.k+=q
                end
            else
                line_box[index_in-1].period[1]=τ_L
                line_box[index_in-1].period[2]=τ_R
                for i in 1:2order+1
                    line_tem=line_box[i]
                    line_tem.index=i
                    line_tem.k+=q
                end
                line_tem=line_box[index_in-1]
                line_tem.k-=q
            end
        end

        return false
    else
        #println("insert_r:"*string(r))
        diagram.order+=1

        if index_out-index_in==2
            line_box[index_in].covered=true
        else
            if !cross_over
                line_box[index_in].covered=false
                line_box[index_out-2].covered=false
            else
                line_box[index_in-1].covered=false
                line_box[index_out-1].covered=false
            end
        end

        if !cross_over
            line_tem=Line(k_in,[τ_L,τ_1], m, μ, index_in, false)
            insert!(line_box, index_in, line_tem)
            line_tem=Line(k_out,[τ_2,τ_R_2], m, μ, index_out, false)
            insert!(line_box, index_out, line_tem)
        else
            if index_out!=index_in
                line_tem=Line(k_out,[τ_2,τ_R_2], m, μ, index_out, false)
                insert!(line_box, index_out, line_tem)
                line_tem=Line(k_in,[τ_L,τ_1], m, μ, index_in, false)
                insert!(line_box, index_in, line_tem)
            else
                line_tem=Line(k_in-q,[τ_L,τ_2], m, μ, index_in-1, false)
                insert!(line_box, index_out-1, line_tem)
                line_tem=Line(k_in-q,[τ_1,τ_R], m, μ, index_in+1, false)
                insert!(line_box, index_in+1, line_tem)
            end
        end

        sign_box=diagram.sign_box
        if index_out-index_in==2
            sign_to_add=[[copy(sign_box[index_in][1]),-1],[-1,1],[1,copy(sign_box[index_in][2])]]
            deleteat!(sign_box, index_in)
            for i in 1:3
                insert!(sign_box, index_in, sign_to_add[4-i])
            end
            
        elseif index_out-index_in>0

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

        elseif index_out-index_in==0
            sign_to_add=[[copy(sign_box[index_in-1][1]),1],[1,-1],[-1,copy(sign_box[index_in-1][2])]]
            deleteat!(sign_box, index_in-1)
            for i in 1:3
                insert!(sign_box, index_in-1, sign_to_add[4-i])
            end

        else
            sign_to_add=[[copy(sign_box[index_out-1][1]),1],[1,copy(sign_box[index_out-1][2])]]
            deleteat!(sign_box, index_out-1)
            for i in 1:2
                insert!(sign_box, index_out-1, sign_to_add[3-i])
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

        if !cross_over
            for arc in diagram.arc_box
                arc_index_in = arc.index_in
                arc_index_out = arc.index_out
                if arc_index_out<=index_in
                    continue
                elseif arc_index_in>=index_out-2
                    arc.index_in+=2
                    arc.index_out+=2
                elseif arc_index_in<index_in && arc_index_out>index_out-2
                    arc.index_out+=2
                elseif arc_index_in>=index_in && arc_index_out<=index_out-2
                    arc.index_in+=1
                    arc.index_out+=1
                elseif index_in<=arc_index_in<index_out-2
                    arc.index_in+=1
                    arc.index_out+=2
                elseif index_in<arc_index_out<=index_out-2
                    arc.index_out+=1
                end
            end

            for arc in diagram.end_arc_box
                arc_index_in = arc.index_in
                arc_index_out = arc.index_out
                if arc_index_in<index_in
                    continue
                elseif arc_index_out>index_out-2
                    arc.index_in+=2
                    arc.index_out+=2
                elseif arc_index_in>=index_out-2 && arc_index_out<=index_in
                    arc.index_in+=2
                elseif arc_index_in<index_out-2 && arc_index_out>index_in
                    arc.index_in+=1
                    arc.index_out+=1
                elseif arc_index_out>index_in && arc_index_in>=index_out-2
                    arc.index_in+=2
                    arc.index_out+=1
                elseif index_in>=arc_index_out && arc_index_in<index_out-2
                    arc.index_in+=1
                end
            end
        else
            for arc in diagram.arc_box
                arc_index_in = arc.index_in
                arc_index_out = arc.index_out
                if arc_index_out<=index_out-1
                    continue
                elseif arc_index_in>=index_in-1
                    arc.index_in+=2
                    arc.index_out+=2
                elseif arc_index_in>=index_out-1 && arc_index_out<=index_in-1
                    arc.index_in+=1
                    arc.index_out+=1
                elseif arc_index_in<index_out-1 && arc_index_out<=index_in-1
                    arc.index_out+=1
                elseif arc_index_in<index_out-1 && arc_index_out>index_in-1
                    arc.index_out+=2
                elseif arc_index_in>=index_out-1 && arc_index_out>index_in-1
                    arc.index_in+=1
                    arc.index_out+=2
                end
            end

            for arc in diagram.end_arc_box
                arc_index_in = arc.index_in
                arc_index_out = arc.index_out
                if arc_index_in<index_out-1
                    continue
                elseif arc_index_out>index_in-1
                    arc.index_in+=2
                    arc.index_out+=2
                elseif arc_index_in>=index_in-1 && arc_index_out<=index_out-1
                    arc.index_in+=2
                elseif arc_index_in<index_in-1 && arc_index_out<=index_out-1
                    arc.index_in+=1
                elseif arc_index_in<index_in-1 && arc_index_out>index_out-1
                    arc.index_in+=1
                    arc.index_out+=1
                elseif arc_index_in>=index_in-1 && arc_index_out>index_out-1
                    arc.index_out+=1
                    arc.index_in+=2
                end
            end
        end

        if !cross_over
            push!(diagram.arc_box,new_arc)
        else
            push!(diagram.end_arc_box,new_arc)
        end


        return true
    end
end

function remove_arc!(diagram::Diagram,order::Int64,m::Float64,μ::Float64,ω::Float64,α_squared::Float64)

    if order-1<0
        return false
    end
    
    arc_box=diagram.arc_box
    end_arc_box = diagram.end_arc_box
    index=rand(1:order)
    τ=0
    if index <= length(arc_box)
        arc=arc_box[index]
        closed_arc = true
    else
        arc=end_arc_box[index-length(arc_box)]
        closed_arc = false
        τ=diagram.τ
    end

    index_in=arc.index_in
    index_out=arc.index_out
    q=arc.q


    line_box=diagram.line_box
    line_in=line_box[index_in]
    line_out=line_box[index_out]


    τ_L=line_in.period[1]

    if index_out-index_in==2
        τ_R=line_out.period[2]
    elseif index_out==index_in
        τ_L=line_box[index_in-1].period[1]
        τ_R=line_box[index_in+1].period[2]
    else
        τ_R=line_box[index_in+1].period[2]
    end

    τ_R_2=line_out.period[2]
    
    τ_1=arc.period[1]
    τ_2=arc.period[2]
    arc_T=abs(τ-abs(τ_2-τ_1))
    total_dis=0
    w_x=1
    w_y=exp(-ω*arc_T)*α_squared/(2*pi)^3/norm(q)^2
    # println("τ"*string(τ))
    if closed_arc
        for i in index_in+1:index_out-1
            line_tem=line_box[i]
            total_dis+=dispersion(line_tem)
            line_tem.index-=1
            line_tem.k+=q
            total_dis-=dispersion(line_tem)
        end
    else
        for i in [collect(1:index_out-1); collect(index_in+1:2order+1)]
            line_tem=line_box[i]
            total_dis+=dispersion(line_tem)
            line_tem.k+=q
            total_dis-=dispersion(line_tem)
        end
        for i in index_out:index_in
            line_tem=line_box[i]
            line_tem.index-=1
        end
        for i in index_in+1:2order+1
            line_tem=line_box[i]
            line_tem.index-=2
        end
            
    end

    τ=diagram.τ
    p_x_y=diagram.p_ins/(2order-1)/(τ_R-τ_L)
    p_x_y*=exp(-norm(q)^2/(2m)*arc_T)/(2pi*m/arc_T)^1.5

    p_y_x=diagram.p_rem/order#*(1-exp(-ω*τ))/ω
    r=((2*pi)^3*p_x_y*norm(arc.q)^2)/(exp(total_dis)*p_y_x*α_squared)
    
    if r<rand()
        if closed_arc
            for i in index_in+1:index_out-1
                line_tem=line_box[i]
                line_tem.index+=1
                line_tem.k-=q
            end
        else
            for i in [collect(1:index_out-1); collect(index_in+1:2order+1)]
                line_tem=line_box[i]
                line_tem.k-=q
            end
            for i in index_out:index_in
                line_tem=line_box[i]
                line_tem.index+=1
            end
            for i in index_in+1:2order+1
                line_tem=line_box[i]
                line_tem.index+=2
            end
        end    

        return false
    else
        #println("delete_r:"*string(r))
        sign_box=diagram.sign_box
        diagram.order-=1
        if closed_arc
            deleteat!(arc_box, index)
        else
            deleteat!(end_arc_box, index-length(arc_box))
        end

        if closed_arc
            if index_out-index_in==2
                line_tem=line_box[index_in+1]
                line_tem.period[1]=τ_L
                line_tem.period[2]=τ_R
                line_tem.covered=false
                sign_to_add=[sign_box[index_in][1],sign_box[index_out][2]]
                deleteat!(sign_box, index_in:index_out)
                insert!(sign_box, index_in, sign_to_add)

            else
                line_tem=line_box[index_in+1]
                line_tem.period[1]=τ_L
                sign_to_add=[sign_box[index_in][1],sign_box[index_in+1][2]]
                deleteat!(sign_box, index_in:index_in+1)
                insert!(sign_box, index_in, sign_to_add)


                line_tem=line_box[index_out-1]
                line_tem.period[2]=τ_R_2
                sign_to_add=[sign_box[index_out-2][1],sign_box[index_out-1][2]]
                deleteat!(sign_box, index_out-2:index_out-1)
                insert!(sign_box, index_out-2, sign_to_add)

            end
        else
            if index_out == index_in
                line_tem=line_box[index_in]
                line_tem.period[1]=τ_L
                line_tem.period[2]=τ_R
                line_tem.covered=false
                sign_to_add=[sign_box[index_out-1][1],sign_box[index_in+1][2]]
                deleteat!(sign_box, index_out-1:index_in+1)
                insert!(sign_box, index_out-1, sign_to_add)

            else
                line_tem=line_box[index_in]
                line_tem.period[1]=τ_L
                line_tem.period[2]=τ_R
                sign_to_add=[sign_box[index_in][1],sign_box[index_in+1][2]]
                deleteat!(sign_box,index_in:index_in+1)
                insert!(sign_box, index_in, sign_to_add)

                line_tem=line_box[index_out]
                line_tem.period[1]=line_box[index_out-1].period[1]
                sign_to_add=[sign_box[index_out-1][1],sign_box[index_out][2]]
                deleteat!(sign_box,index_out-1:index_out)
                insert!(sign_box, index_out-1, sign_to_add)

            end
        end

        if closed_arc
            deleteat!(line_box, [index_in,index_out])
        else
            deleteat!(line_box, [index_out-1, index_in+1])
        end

        for i in 1:length(line_box)
            line_box[i].index=i
        end

        if closed_arc
            for arc in diagram.arc_box
                arc_index_in = arc.index_in
                arc_index_out = arc.index_out
                if arc_index_out<=index_in
                    continue
                elseif arc_index_in>=index_out
                    arc.index_in-=2
                    arc.index_out-=2
                elseif arc_index_in<index_in && arc_index_out>index_out
                    arc.index_out-=2
                    if arc.index_out-arc.index_in == 2
                        line_box[arc.index_in+1].covered=true
                    end
                elseif arc_index_in>index_in && arc_index_out<index_out
                    arc.index_in-=1
                    arc.index_out-=1
                elseif index_in<arc_index_in<index_out
                    arc.index_in-=1
                    arc.index_out-=2
                    if arc.index_out-arc.index_in == 2
                        line_box[arc.index_in+1].covered=true
                    end
                elseif index_in<arc_index_out<index_out
                    arc.index_out-=1
                    if arc.index_out-arc.index_in == 2
                        line_box[arc.index_in+1].covered=true
                    end
                end
            end

            for arc in diagram.end_arc_box
                arc_index_in = arc.index_in
                arc_index_out = arc.index_out
                if arc_index_in < index_in
                    continue
                elseif arc_index_out > index_out
                    arc.index_in-=2
                    arc.index_out-=2
                elseif arc_index_out <= index_in && arc_index_in >= index_out
                    arc.index_in-=2
                elseif arc_index_out > index_in && arc_index_in < index_out
                    arc.index_in-=1
                    arc.index_out-=1
                elseif index_in < arc_index_out < index_out
                    arc.index_out-=1
                    arc.index_in-=2
                elseif index_in < arc_index_in < index_out
                    arc.index_in-=1
                end
            end

        else
            for arc in diagram.arc_box
                arc_index_in = arc.index_in
                arc_index_out = arc.index_out
                if arc_index_out < index_out
                    continue
                elseif arc_index_in > index_in
                    arc.index_in-=2
                    arc.index_out-=2
                elseif arc_index_in < index_out && arc_index_out > index_in
                    arc.index_out-=2
                    if arc.index_out-arc.index_in == 2
                        line_box[arc.index_in+1].covered=true
                    end
                elseif arc_index_in >= index_out && arc_index_out <= index_in
                    arc.index_in-=1
                    arc.index_out-=1
                elseif index_out <= arc_index_in < index_in
                    arc.index_in-=1
                    arc.index_out-=2
                    if arc.index_out-arc.index_in == 2
                        line_box[arc.index_in+1].covered=true
                    end
                elseif index_out < arc_index_out <= index_in
                    arc.index_out-=1
                    if arc.index_out-arc.index_in == 2
                        line_box[arc.index_in+1].covered=true
                    end
                end
            end
            
            for arc in diagram.end_arc_box
                arc_index_in = arc.index_in
                arc_index_out = arc.index_out
                if arc_index_in < index_out
                    continue
                elseif arc_index_out > index_in
                    arc.index_in-=2
                    arc.index_out-=2
                elseif arc_index_out < index_out && arc_index_in > index_in
                    arc.index_in-=2
                elseif arc_index_out > index_out && arc_index_in < index_in
                    arc.index_in-=1
                    arc.index_out-=1
                elseif index_out < arc_index_out <= index_in
                    arc.index_in-=2
                    arc.index_out-=1
                elseif index_out <= arc_index_in < index_in
                    arc.index_in-=1
                end
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
    end_arc_box=diagram.end_arc_box
    left_check=false
    right_check=false
    left_open=false
    right_open=false
    left_index=0
    right_index=0

    for i in 1:length(arc_box)
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

    for i in 1:length(end_arc_box)
        arc=end_arc_box[i]
        if !left_check
            if arc_judge(arc,sign[1],false,line_index)
                left_index=i
                left_check=true
                left_open=true
            end
        end

        if !right_check
            if arc_judge(arc,sign[2],true,line_index)
                right_index=i
                right_check=true
                right_open=true
            end
        end

        if right_check && left_check
            break
        end
    end

    if right_open && left_open && right_index == left_index
        return false
    end

    # println("swap_index is:",line_index)
    if left_open
        arc_l=end_arc_box[left_index]
    else
        arc_l=arc_box[left_index]
    end
    
    if right_open
        arc_r=end_arc_box[right_index]
    else
        arc_r=arc_box[right_index]
    end

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

        if left_open
            deleteat!(end_arc_box, left_index)
            insert!(end_arc_box, left_index, new_arc_l)
        else
            deleteat!(arc_box, left_index)
            insert!(arc_box, left_index, new_arc_l)
        end

        
        if right_open
            deleteat!(end_arc_box, right_index)
            insert!(end_arc_box, right_index, new_arc_r)
        else
            deleteat!(arc_box, right_index)
            insert!(arc_box, right_index, new_arc_r)
        end

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

    line_box=diagram.line_box
    order=diagram.order
    line_end=line_box[2*order+1]
    p=line_end.k
    μ=diagram.μ
    m=diagram.mass
    ω=diagram.ω
    n_phonon=length(diagram.end_arc_box)
    dispersion=norm(p)^2/(2m)-μ+n_phonon*ω

    τ_new=line_end.period[1]-log(rand())/(dispersion)

    if τ_new>diagram.max_τ
        return false
    end

    line_end.period[2]=τ_new
    line_box[2*order+1]=line_end
    diagram.τ=τ_new

    return true

end


function scale!(diagram::Diagram)

    line_box=diagram.line_box
    order=diagram.order
    arc_box=diagram.arc_box
    end_arc_box=diagram.end_arc_box
    μ=diagram.μ
    m=diagram.mass
    ω=diagram.ω
    τ=diagram.τ
    record_τ=τ=diagram.record_τ
    total_dis=0

    for line in line_box
        total_dis+=dispersion(line)
    end

    for arc in arc_box
        period=arc.period
        total_dis+=-ω*(period[2]-period[1])
    end

    for arc in end_arc_box
        period=arc.period
        total_dis+=-ω*(τ+period[2]-period[1])
    end

    coef=-τ/total_dis
    #τ_new=-log(rand())/(guess_E-μ)
    τ_new=rand(Gamma(2*order+1,coef))
    #println(coef)
    
    if τ_new>diagram.max_τ
        return false
    end

    # #τ_new=(1+λ)*record_τ
    # λ=τ_new/record_τ-1
    # total_dis=total_dis*λ
    # total_dis+=2*order*log(1+λ)
    # println(2*order*log(1+λ))
    
    # total_dis+=(guess_E-μ)*(λ*record_τ)
    # # y_x=exp((guess_E-μ)*(λ*record_τ))
    # # x_y=exp(-(guess_E-μ)*τ_new)
    
    # r=exp(total_dis)

    # if r<rand()
    #     return false
    # else
    diagram.record_τ=τ_new
    return true

end
