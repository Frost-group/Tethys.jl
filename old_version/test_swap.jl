using SpecialFunctions
using Roots

function swap_arc!(diagram::Diagram)

    order=diagram.order
    m=diagram.mass
    μ=diagram.μ
    ω=diagram.ω
    τ=diagram.τ

    if order<2
        return false
    end
    # println(diagram.line_box)
    line_box=diagram.line_box
    line_index=sort(sample(1:2*order,2,replace = false))
    # println(line_index)
    index_l=line_index[1]+1
    index_r=line_index[2]
    # chosen_line=line_box[line_index]
    chose_l=line_box[index_l]
    chose_r=line_box[index_r]
    
    if index_l==index_r && chose_l.covered
        return false
    end

    sign_l=Int(diagram.sign_box[index_l][1])
    sign_r=Int(diagram.sign_box[index_r][2])
    arc_box=diagram.arc_box
    end_arc_box=diagram.end_arc_box
    arc_box_length = length(arc_box)
    left_check=false
    right_check=false
    left_open=false
    right_open=false
    left_index=0
    right_index=0

    for i in 1:arc_box_length
        arc=arc_box[i]
        if !left_check
            if arc_judge(arc,sign_l,false,index_l)
                left_index=i
                left_check=true
            end
        end

        if !right_check
            if arc_judge(arc,sign_r,true,index_r)
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
            if arc_judge(arc,sign_l,false,index_l)
                left_index=i
                left_check=true
                left_open=true
            end
        end

        if !right_check
            if arc_judge(arc,sign_r,true,index_r)
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

    if sign_l == 1
        new_arc_l=Arc(q1,[arc_l.period[1],chose_r.period[2]],ω,arc_l.index_in,index_r+1)
    else
        new_arc_l=Arc(q1,[chose_r.period[2],arc_l.period[2]],ω,index_r,arc_l.index_out)
    end

    if sign_r == 1
        new_arc_r=Arc(q2,[arc_r.period[1],chose_l.period[1]],ω,arc_r.index_in,index_l)
    else
        new_arc_r=Arc(q2,[chose_l.period[1],arc_r.period[2]],ω,index_l-1,arc_r.index_out)
    end
    
    total_dis=0

    for i in index_l:index_r
        line_tem=line_box[i]
        total_dis-=dispersion(line_tem, m, μ)*τ
        line_tem.k+=-sign_l*q1+sign_r*q2
        total_dis+=dispersion(line_tem, m, μ)*τ
    end
    # new_line=Line(chosen_line.k-sign[1]*q1+sign[2]*q2 ,chosen_line.period, line_index,false)

    # total_dis=dispersion(new_line, m, μ)-dispersion(chosen_line, m, μ)
    total_dis+=(arc_dispersion(new_arc_r,ω)+arc_dispersion(new_arc_l,ω)-arc_dispersion(arc_r,ω)-arc_dispersion(arc_l,ω))*τ
    r=exp(total_dis)
    # println(w_y/w_x)
    # r=ratio*(1+log(ratio)/diagram.dispersion)^(2order+1)
    
    # println("          ")
    if r<rand()
        # println(line_box)
        # println(index_l)
        # println(index_r)
        for i in index_l:index_r
            line_tem=line_box[i]
            line_tem.k-=-sign_l*q1+sign_r*q2
        end
        println(line_box)
        return false
        diagram.dispersion+=total_dis
        sign_box=diagram.sign_box

        if index_l!=index_r
            sign_box[index_l]=[sign_r,sign_box[index_l][2]]
            sign_box[index_r]=[sign_box[index_r][1],sign_l]
        else
            sign_box[index_l]=[sign_r,sign_l]
        end
        sign_box[index_l-1]=[sign_box[index_l-1][1],sign_r]
        sign_box[index_r+1]=[sign_l,sign_box[index_r+1][2]]

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



function resample_arc!(diagram::Diagram,order::Int64,m::Int64,μ::Float64,ω::Int64,α_squared::Float64)

    if order-1<0
        return false
    end
    
    arc_box=diagram.arc_box
    end_arc_box = diagram.end_arc_box
    arc_box_length = length(arc_box)
    index=rand(1:order)
    offset_τ=0
    τ=diagram.τ
    println(index)
    if index <= arc_box_length
        arc=arc_box[index]
        closed_arc = true
        
    else
        arc=end_arc_box[index-arc_box_length]
        closed_arc = false
        offset_τ=diagram.τ
    end

    index_in=arc.index_in
    index_out_r=arc.index_out
    q=arc.q


    line_box=diagram.line_box
    line_in=line_box[index_in]
    line_out=line_box[index_out_r]

    τ_1=arc.period[1]*τ
    τ_2=arc.period[2]*τ
    arc_T=abs(offset_τ-abs(τ_2-τ_1))
    total_dis=0
    # w_x=1
    w_y=exp(-ω*arc_T)*α_squared/(2*pi)^3/norm(q)^2

    open_arc_range = [collect(1:index_out_r-1); collect(index_in+1:2order+1)]

    if closed_arc
        for i in index_in+1:index_out_r-1
            line_tem=line_box[i]
            total_dis-=dispersion(line_tem, m, μ)*τ
            line_tem.index-=1
            line_tem.k+=q
            total_dis+=dispersion(line_tem, m, μ)*τ
        end
    else
        for i in open_arc_range
            line_tem=line_box[i]
            total_dis-=dispersion(line_tem, m, μ)*τ
            line_tem.k+=q
            total_dis+=dispersion(line_tem, m, μ)*τ
        end        
    end

    τ_2_in=τ_1-log(rand())/ω
    arc_T_in=τ_2_in-τ_1
    if arc_T_in > τ
        println("overlap")
        return false
    end

    phi = rand(Uniform(0,pi*2))
    costheta = rand(Uniform(-1,1))
    theta = acos(costheta)
    x = sin(theta)*cos(phi)
    y = sin(theta)*sin(phi)
    z = cos(theta)
    q_in = abs(rand(Normal(0,sqrt(m/arc_T_in)))).*[x,y,z]

    w_x=1.0
    w_y=1.0
    total_dis=0.0
    τ_R_2=0.0
    index_out_in=0
    k_out=0

    if τ_2_in<τ
        #not set covered yet
        for i in index_in:2order+1
            line_tem=line_box[i]
            if line_tem.period[2]<τ_2_in/τ
                total_dis-=dispersion(line_tem, m, μ)*τ
                line_tem.k-=q_in
                line_tem.index+=1
                total_dis+=dispersion(line_tem, m, μ)*τ
                continue
            else
                k_out=deepcopy(line_tem.k)
                τ_R_2=deepcopy(line_tem.period[2])
                line_tem.period[2]=τ_2_in/τ
                total_dis-=dispersion(line_tem, m, μ)*τ
                line_tem.k-=q_in
                line_tem.index+=1
                total_dis+=dispersion(line_tem, m, μ)*τ
                index_out_in=i+1
                break
            end
        end
    else
        cross_over=true
        τ_2_in=τ_2_in-τ
        for i in index_in:2order+1
            line_tem=line_box[i]
            total_dis-=dispersion(line_tem, m, μ)*τ
            line_tem.k-=q_in
            line_tem.index+=1
            total_dis+=dispersion(line_tem, m, μ)*τ
        end

        for i in 1:index_in
            line_tem=line_box[i]
            if line_tem.period[2]<τ_2_in/τ
                total_dis-=dispersion(line_tem, m, μ)*τ
                line_tem.k-=q_in
                total_dis+=dispersion(line_tem, m, μ)*τ
                continue
            else
                k_out=deepcopy(line_tem.k)
                τ_R_2=deepcopy(line_tem.period[2])
                line_tem.period[2]=τ_2_in/τ
                total_dis-=dispersion(line_tem, m, μ)*τ
                line_tem.k-=q_in
                total_dis+=dispersion(line_tem, m, μ)*τ
                index_out_in=i+1
                break
            end
        end
        index_in=index_in+1
    end


    new_arc=Arc(q,[τ_1,τ_2_in]/τ,ω,index_in,index_out_in)

    p_x_y=exp(-norm(q)^2/(2m)*arc_T)/(2pi*m/arc_T)^0.5*exp(-ω*arc_T)
    p_x_y/=exp(-norm(q_in)^2/(2m)*arc_T_in)/(2pi*m/arc_T_in)^0.5*exp(-ω*arc_T_in)
    r=exp(total_dis)*p_x_y

    if r<rand()

        if !cross_over
            line_box[index_out_in-1].period[2]=τ_R_2
            for i in index_in:index_out_in-2
                line_tem=line_box[i]
                line_tem.index=i
                line_tem.k+=q_in
            end
        else
            line_box[index_out_in-1].period[2]=τ_R_2
            for i in index_in-1:2order+1
                line_tem=line_box[i]
                line_tem.index=i
                line_tem.k+=q_in
            end

            for i in 1:index_out_in-1
                line_tem=line_box[i]
                line_tem.index=i
                line_tem.k+=q_in
            end
        end

        if closed_arc
            for i in index_in+1:index_out_r-1
                line_tem=line_box[i]
                line_tem.index+=1
                line_tem.k-=q
            end
        else
            for i in open_arc_range
                line_tem=line_box[i]
                line_tem.k-=q
            end
        end    

        return false

    else

        sign_box=diagram.sign_box
        diagram.dispersion+=total_dis
        diagram.dispersion+=(arc_T-arc_T_in)*ω 

        if closed_arc
            deleteat!(arc_box, index)
        else
            deleteat!(end_arc_box, index-arc_box_length)
        end

        if index_out_r == index_out_in
            line_tem=line_box[index_out_r]
            line_tem.period[1]=τ_2_in/τ
        elseif index_out_r == index_out_in-1
            line_tem_1=line_box[index_out_r-1]
            line_tem_2=line_box[index_out_r]
            τ_c=line_tem_2.period[2]
            line_tem_1.period[2]=τ_c
            deleteat!(line_box, index_out_r)
            line_tem_3=Line(k_out,[τ_c,τ_R_2], index_out_r, false)
            insert!(line_box, index_out_r, line_tem_3)
            new_arc.index_out=index_out_r
        else
            line_tem_1=line_box[index_out_r-1]
            line_tem_2=line_box[index_out_r]
            τ_c=line_tem_2.period[2]
            line_tem_1.period[2]=τ_c

            line_tem_3=Line(k_out,[τ_2_in/τ,τ_R_2], index_out_in, false)
            # insert!(line_box, index_out_in, line_tem_3)

            deleteat!(line_box, index_out_r)
            sign_to_add=[sign_box[index_out_r-1][1],sign_box[index_out_r][2]]
            deleteat!(sign_box, index_out_r-1:index_out_r)
            insert!(sign_box, index_out_r-1, sign_to_add)

            if index_out_in>index_out_r
                index_out_in-=1
                new_arc.index_out=index_out_in

                sign_to_add=[[sign_box[index_out_in-1][1],1],[1,sign_box[index_out_in][2]]]
                deleteat!(sign_box,index_out_in-1:index_out_in)

                for i in 1:2
                    insert!(sign_box, index_out_in-1, sign_to_add[3-i])
                end

            else
                new_arc.index_in+=1

                sign_to_add=[[sign_box[index_out_in-1][1],1],[1,sign_box[index_out_in][2]]]
                deleteat!(sign_box,index_out_in-1:index_out_in)

                for i in 1:2
                    insert!(sign_box, index_out_in-1, sign_to_add[3-i])
                end
            end

        end
        #delete

        index_in=arc.index_in
        index_out=arc.index_out

        if closed_arc
            for arc in diagram.arc_box
                # total_dispersion+=arc_dispersion(arc,ω)
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
                # total_dispersion+=end_arc_dispersion(arc,ω, τ)
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
                # total_dispersion+=arc_dispersion(arc,ω)
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
                # total_dispersion+=end_arc_dispersion(arc,ω, τ)
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
            diagram.component-=1
        end

        #insert
        index_in=new_arc.index_in
        index_out=new_arc.index_out
        if !cross_over
            for arc in diagram.arc_box
                # total_dispersion+=arc_dispersion(arc,ω)
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
                # total_dispersion+=end_arc_dispersion(arc,ω, τ)
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
                # total_dispersion+=arc_dispersion(arc,ω)
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
                # total_dispersion+=end_arc_dispersion(arc,ω, τ)
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
            diagram.component+=1
        end

        return true
    end
end

function shift_2!(diagram::Diagram,sign::Int64,line_index::Int64,arc_index::Int64,closed_arc::Bool,m::Int64,μ::Float64)

    arc_box=diagram.arc_box
    end_arc_box = diagram.end_arc_box
    line_box=diagram.line_box
    ω=diagram.ω
    τ=diagram.τ
    old_length=0

    if closed_arc
        arc=arc_box[arc_index]
        old_length=(arc.period[2]-arc.period[1])*τ
    else
        arc=end_arc_box[arc_index]
        old_length=τ+(arc.period[2]-arc.period[1])*τ
    end

    if sign==-1
        old_τ=arc.period[1] 
        unchanged_τ=arc.period[2] 
    else
        old_τ=arc.period[2]
        unchanged_τ=arc.period[1] 
    end

    line=line_box[line_index]
    line_next=line_box[line_index+1]  

    τ_a=line.period[1]*τ
    τ_c=line_next.period[2]*τ
    # τ_2=τ_1-log(rand())/ω
    # e=(norm(line.k)^2-norm(line_next.k)^2)/(2m)+sign*ω
    # ω/(1-exp(-ω*(τ_c-τ_a)))
    if sign==1
        τ_b=unchanged_τ*τ-log(rand())/ω
    else
        τ_b=unchanged_τ*τ+log(rand())/ω
    end

    if !closed_arc
        τ_b=τ_b-τ
    end

    if τ_b<τ_a || τ_b>τ_c
        # println("reject")
        return false
    end

    # if isnan(τ_b)
    #     return false
    # end

    total_dis=0
    total_dis-=dispersion(line, m, μ)*τ+dispersion(line_next, m, μ)*τ
    total_dis-=arc_dispersion(arc, ω)*τ#(-1)^(!closed_arc)*

    line.period[2]=τ_b/τ
    line_next.period[1]=τ_b/τ
    # total_dis=0

    if sign==-1
        arc.period[1]=τ_b/τ
    else
        arc.period[2]=τ_b/τ
    end
    total_dis+=dispersion(line, m, μ)*τ+dispersion(line_next, m, μ)*τ
    total_dis+=arc_dispersion(arc, ω)*τ#(-1)^(!closed_arc)*
    r=exp(total_dis)

    if closed_arc
        arc=arc_box[arc_index]
        new_length=(arc.period[2]-arc.period[1])*τ
    else
        arc=end_arc_box[arc_index]
        new_length=τ+(arc.period[2]-arc.period[1])*τ
    end

    r*=exp(-ω*(old_length))/exp(-ω*(new_length))

    # if sign==1
    #     r*=exp(-ω*(old_length))/exp(-ω*(τ_b-τ_a))
    # else
    #     r*=exp(-ω*(τ_c-old_τ*τ))/exp(-ω*(τ_c-τ_b))
    # end
    
    if r<rand()
        if sign==-1
            arc.period[1]=old_τ
        else
            arc.period[2]=old_τ
        end

        line.period[2]=old_τ
        line_next.period[1]=old_τ
        return false
    else
        diagram.dispersion+=total_dis
        return true
    end

end

function shift_3!(diagram::Diagram,sign::Int64,line_index::Int64,arc_index::Int64,closed_arc::Bool,m::Int64,μ::Float64,order::Int64)

    arc_box=diagram.arc_box
    end_arc_box = diagram.end_arc_box
    line_box=diagram.line_box
    ω=diagram.ω
    τ=diagram.τ
    old_length=0
    offset_τ=0


    if closed_arc
        arc=arc_box[arc_index]
        old_length=(arc.period[2]-arc.period[1])*τ
    else
        arc=end_arc_box[arc_index]
        old_length=τ+(arc.period[2]-arc.period[1])*τ
        offset_τ=diagram.τ
    end

    if sign==-1
        old_τ=arc.period[1] 
        unchanged_τ=arc.period[2] 
    else
        old_τ=arc.period[2]
        unchanged_τ=arc.period[1] 
    end

    line=line_box[line_index]
    line_next=line_box[line_index+1]  

    τ_a=line.period[1]*τ
    τ_c=line_next.period[2]*τ
    # τ_2=τ_1-log(rand())/ω
    # e=(norm(line.k)^2-norm(line_next.k)^2)/(2m)+sign*ω
    # ω/(1-exp(-ω*(τ_c-τ_a)))
    if sign==1
        τ_b=unchanged_τ*τ-log(rand())/ω
    else
        τ_b=unchanged_τ*τ+log(rand())/ω
    end

    if !closed_arc
        τ_b=τ_b-τ
    end

    if τ_b<τ_a || τ_b>τ_c
        # println("reject")
        return false
    end

    # if isnan(τ_b)
    #     return false
    # end

    total_dis=0
    total_dis-=dispersion(line, m, μ)*τ+dispersion(line_next, m, μ)*τ
    total_dis-=arc_dispersion(arc, ω)*τ#(-1)^(!closed_arc)*

    line.period[2]=τ_b/τ
    line_next.period[1]=τ_b/τ
    # total_dis=0

    if sign==-1
        arc.period[1]=τ_b/τ
    else
        arc.period[2]=τ_b/τ
    end
    total_dis+=dispersion(line, m, μ)*τ+dispersion(line_next, m, μ)*τ
    total_dis+=arc_dispersion(arc, ω)*τ#(-1)^(!closed_arc)*
    r=1

    if closed_arc
        arc=arc_box[arc_index]
        new_length=(arc.period[2]-arc.period[1])*τ
    else
        arc=end_arc_box[arc_index]
        new_length=τ+(arc.period[2]-arc.period[1])*τ
    end

    r*=exp(-ω*(old_length))/exp(-ω*(new_length))

    index_in=arc.index_in
    index_out=arc.index_out
    q=arc.q

    τ_1=arc.period[1]*τ
    τ_2=arc.period[2]*τ
    arc_T=abs(offset_τ-abs(τ_2-τ_1))

    phi = rand(Uniform(0,pi*2))
    costheta = rand(Uniform(-1,1))
    theta = acos(costheta)
    x = sin(theta)*cos(phi)
    y = sin(theta)*sin(phi)
    z = cos(theta)
    p = abs(rand(Normal(0,sqrt(m/arc_T)))).*[x,y,z]

    open_arc_range = [collect(1:index_out-1); collect(index_in+1:2order+1)]

    if closed_arc
        for i in index_in+1:index_out-1
            line_tem=line_box[i]
            total_dis-=dispersion(line_tem, m, μ)*τ
            line_tem.k+=q
            line_tem.k-=p
            total_dis+=dispersion(line_tem, m, μ)*τ
        end
    else
        for i in open_arc_range
            line_tem=line_box[i]
            total_dis-=dispersion(line_tem, m, μ)*τ
            line_tem.k+=q
            line_tem.k-=p
            total_dis+=dispersion(line_tem, m, μ)*τ
        end       
    end

    r=exp(total_dis)
    r*=exp(-(norm(q)^2-norm(p)^2)/(2m)*arc_T)
    
    if r<rand()

        if closed_arc
            for i in index_in+1:index_out-1
                line_tem=line_box[i]
                line_tem.k+=p
                line_tem.k-=q
            end
        else
            for i in open_arc_range
                line_tem=line_box[i]
                line_tem.k+=p
                line_tem.k-=q
            end       
        end

        if sign==-1
            arc.period[1]=old_τ
        else
            arc.period[2]=old_τ
        end

        line.period[2]=old_τ
        line_next.period[1]=old_τ
        return false
    else
        arc.q=p
        diagram.dispersion+=total_dis
        return true
    end

end

function update_arcp!(diagram::Diagram,order::Int64,m::Int64,μ::Float64)
    
    if order-1<0
        return false
    end

    arc_box=diagram.arc_box
    end_arc_box = diagram.end_arc_box
    arc_box_length = length(arc_box)
    index=rand(1:order)
    offset_τ=0
    τ=diagram.τ
    if index <= arc_box_length
        arc=arc_box[index]
        closed_arc = true
        
    else
        arc=end_arc_box[index-arc_box_length]
        closed_arc = false
        offset_τ=diagram.τ
    end

    index_in=arc.index_in
    index_out=arc.index_out
    q=arc.q

    line_box=diagram.line_box

    τ_1=arc.period[1]*τ
    τ_2=arc.period[2]*τ
    arc_T=abs(offset_τ-abs(τ_2-τ_1))

    phi = rand(Uniform(0,pi*2))
    costheta = rand(Uniform(-1,1))
    theta = acos(costheta)
    x = sin(theta)*cos(phi)
    y = sin(theta)*sin(phi)
    z = cos(theta)
    p = abs(rand(Normal(0,sqrt(m/arc_T)))).*[x,y,z]

    total_dis=0

    open_arc_range = [collect(1:index_out-1); collect(index_in+1:2order+1)]

    if closed_arc
        for i in index_in+1:index_out-1
            line_tem=line_box[i]
            total_dis+=dispersion(line_tem, m, μ)*τ
            line_tem.k+=q
            line_tem.k-=p
            total_dis-=dispersion(line_tem, m, μ)*τ
        end
    else
        for i in open_arc_range
            line_tem=line_box[i]
            total_dis+=dispersion(line_tem, m, μ)*τ
            line_tem.k+=q
            line_tem.k-=p
            total_dis-=dispersion(line_tem, m, μ)*τ
        end       
    end

    r=exp(-total_dis)
    r*=exp(-(norm(q)^2-norm(p)^2)/(2m)*arc_T)

    if r<rand()
        if closed_arc
            for i in index_in+1:index_out-1
                line_tem=line_box[i]
                line_tem.k+=p
                line_tem.k-=q
            end
        else
            for i in open_arc_range
                line_tem=line_box[i]
                line_tem.k+=p
                line_tem.k-=q
            end       
        end
        return false
    else
        diagram.dispersion-=total_dis
        arc.q=p
        return true
    end
end

function update_arcp!(diagram::Diagram,order::Int64,index::Int64,closed_arc::Bool,m::Int64,μ::Float64)

    arc_box=diagram.arc_box
    end_arc_box = diagram.end_arc_box
    arc_box_length = length(arc_box)
    offset_τ=0
    τ=diagram.τ

    if closed_arc
        arc=arc_box[index]
    else
        arc=end_arc_box[index]
        offset_τ=diagram.τ
    end

    index_in=arc.index_in
    index_out=arc.index_out
    q=arc.q

    line_box=diagram.line_box

    τ_1=arc.period[1]*τ
    τ_2=arc.period[2]*τ
    arc_T=abs(offset_τ-abs(τ_2-τ_1))

    phi = rand(Uniform(0,pi*2))
    costheta = rand(Uniform(-1,1))
    theta = acos(costheta)
    x = sin(theta)*cos(phi)
    y = sin(theta)*sin(phi)
    z = cos(theta)
    p = abs(rand(Normal(0,sqrt(m/arc_T)))).*[x,y,z]

    total_dis=0

    open_arc_range = [collect(1:index_out-1); collect(index_in+1:2order+1)]

    if closed_arc
        for i in index_in+1:index_out-1
            line_tem=line_box[i]
            total_dis+=dispersion(line_tem, m, μ)*τ
            line_tem.k+=q
            line_tem.k-=p
            total_dis-=dispersion(line_tem, m, μ)*τ
        end
    else
        for i in open_arc_range
            line_tem=line_box[i]
            total_dis+=dispersion(line_tem, m, μ)*τ
            line_tem.k+=q
            line_tem.k-=p
            total_dis-=dispersion(line_tem, m, μ)*τ
        end       
    end

    r=exp(-total_dis)
    r*=exp(-(norm(q)^2-norm(p)^2)/(2m)*arc_T)

    if r<rand()
        if closed_arc
            for i in index_in+1:index_out-1
                line_tem=line_box[i]
                line_tem.k+=p
                line_tem.k-=q
            end
        else
            for i in open_arc_range
                line_tem=line_box[i]
                line_tem.k+=p
                line_tem.k-=q
            end       
        end
        return false
    else
        diagram.dispersion-=total_dis
        arc.q=p
        return true
    end
end

function shift!(diagram::Diagram)

    arc_box=diagram.arc_box
    end_arc_box = diagram.end_arc_box
    line_box=diagram.line_box
    ω=diagram.ω
    τ=diagram.τ

    if closed_arc
        arc=arc_box[arc_index]
    else
        arc=end_arc_box[arc_index]
    end

    line=line_box[line_index]
    line_next=line_box[line_index+1]  
    diagram.dispersion-=dispersion(line, m, μ)*τ+dispersion(line_next, m, μ)*τ
    diagram.dispersion-=arc_dispersion(arc, ω)*τ#(-1)^(!closed_arc)*
    τ_a=line.period[1]*τ
    τ_c=line_next.period[2]*τ
    e=(norm(line.k)^2-norm(line_next.k)^2)/(2m)+sign*ω
    if e>0
        τ_b=τ_a-log(1-rand()*(1-exp(-e*(τ_c-τ_a))))/e
    else
        τ_b=τ_c-log(1-rand()*(1-exp(e*(τ_c-τ_a))))/e
    end

    line.period[2]=τ_b/τ
    line_next.period[1]=τ_b/τ
    total_dis=0

    if sign==-1
        arc.period[1]=τ_b/τ
    else
        arc.period[2]=τ_b/τ
    end
    diagram.dispersion+=dispersion(line, m, μ)*τ+dispersion(line_next, m, μ)*τ
    diagram.dispersion+=arc_dispersion(arc, ω)*τ#(-1)^(!closed_arc)*
end

function shift!(diagram::Diagram,sign::Int64,line_index::Int64,arc_index::Int64,closed_arc::Bool,m::Int64,μ::Float64)

    arc_box=diagram.arc_box
    end_arc_box = diagram.end_arc_box
    line_box=diagram.line_box
    ω=diagram.ω
    τ=diagram.τ

    if closed_arc
        arc=arc_box[arc_index]
    else
        arc=end_arc_box[arc_index]
    end

    line=line_box[line_index]
    line_next=line_box[line_index+1]  

    τ_a=line.period[1]*τ
    τ_c=line_next.period[2]*τ
    e=(norm(line.k)^2-norm(line_next.k)^2)/(2m)+sign*ω
    # println("e")
    # println(e)
    # println(line.k)
    # println(line_next.k)
    if e>0
        τ_b=τ_a-log(1-rand()*(1-exp(-e*(τ_c-τ_a))))/e
    else
        τ_b=τ_c-log(1-rand()*(1-exp(e*(τ_c-τ_a))))/e
    end
    if isnan(τ_b)
        return false
    end
    diagram.dispersion-=dispersion(line, m, μ)*τ+dispersion(line_next, m, μ)*τ
    diagram.dispersion-=arc_dispersion(arc, ω)*τ#(-1)^(!closed_arc)*

    line.period[2]=τ_b/τ
    line_next.period[1]=τ_b/τ
    total_dis=0

    if sign==-1
        arc.period[1]=τ_b/τ
    else
        arc.period[2]=τ_b/τ
    end
    diagram.dispersion+=dispersion(line, m, μ)*τ+dispersion(line_next, m, μ)*τ
    diagram.dispersion+=arc_dispersion(arc, ω)*τ#(-1)^(!closed_arc)*

end

function swap_arc!(diagram::Diagram)

    order=diagram.order
    m=diagram.mass
    μ=diagram.μ
    ω=diagram.ω
    τ=diagram.τ

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
    arc_box_length = length(arc_box)
    left_check=false
    right_check=false
    left_open=false
    right_open=false
    left_index=0
    right_index=0

    for i in 1:arc_box_length
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
    # println("he")
    if right_open && left_open && right_index == left_index
        println("overlap")
        return false
    end
    # println(right_check,left_check)
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

    new_line=Line(chosen_line.k-sign[1]*q1+sign[2]*q2 ,chosen_line.period, line_index,false)
    # w_x=green_zero(chosen_line, m, μ)*phonon_propagator(arc_l)*phonon_propagator(arc_r)
    # w_y=green_zero(new_line, m, μ)*phonon_propagator(new_arc_l)*phonon_propagator(new_arc_r)
    # total_dis=dispersion(new_line, m, μ)-dispersion(chosen_line, m, μ)
    # total_dis+=(arc_dispersion(new_arc_r,ω)-arc_dispersion(arc_r,ω))+(arc_dispersion(new_arc_l,ω)-arc_dispersion(arc_l,ω))#(-1)^right_open*#(-1)^left_open*
    dis=(norm(chosen_line.k-sign[1]*q1+sign[2]*q2)^2-norm(chosen_line.k)^2)/(2m)
    dis+=(sign[1]-sign[2])*ω
    dis*=-(chosen_line.period[2]-chosen_line.period[1])*τ
    # p_x_y=insert_prob(new_arc_r,right_open,ω,m,τ)*insert_prob(new_arc_l,left_open,ω,m,τ)
    # p_y_x=insert_prob(arc_r,right_open,ω,m,τ)*insert_prob(arc_l,left_open,ω,m,τ)
    r=exp(dis)#*p_x_y/p_y_x#total_dis*diagram.τ
    # println(r)
    # p_x_y=ω*exp(-ω*(arc_T))*exp(-norm(q)^2/(2m)*arc_T)/(2pi*m/arc_T)^0.5#/norm(q)^2*1.0/(2pi)

    # println(w_y/w_x)
    # r=ratio*(1+log(ratio)/diagram.dispersion)^(2order+1)

    if r<rand()
        return false
    else
        # println(total_dis*diagram.τ)
        diagram.dispersion+=dis#total_dis*diagram.τ#log(r)
        deleteat!(line_box, line_index)
        insert!(line_box, line_index, new_line)
        # diagram.total_dispersion+=dispersion(new_line, m, μ)
        sign_box=diagram.sign_box
        deleteat!(sign_box, line_index)
        insert!(sign_box, line_index, [sign[2],sign[1]])
        sign_box[line_index-1]=[sign_box[line_index-1][1],sign[2]]
        sign_box[line_index+1]=[sign[1],sign_box[line_index+1][2]]

        if left_open
            # diagram.total_dispersion-=end_arc_dispersion(end_arc_box[left_index],ω, diagram.τ)
            deleteat!(end_arc_box, left_index)
            insert!(end_arc_box, left_index, new_arc_l)
            # diagram.total_dispersion+=end_arc_dispersion(new_arc_l,ω, diagram.τ)
        else
            # diagram.total_dispersion-=arc_dispersion(arc_box[left_index],ω)
            deleteat!(arc_box, left_index)
            insert!(arc_box, left_index, new_arc_l)
            # diagram.total_dispersion+=arc_dispersion(new_arc_l,ω)
        end
        if right_open
            # diagram.total_dispersion-=end_arc_dispersion(end_arc_box[right_index],ω, diagram.τ)
            deleteat!(end_arc_box, right_index)
            insert!(end_arc_box, right_index, new_arc_r)
            # diagram.total_dispersion+=end_arc_dispersion(new_arc_r,ω, diagram.τ)
        else
            # diagram.total_dispersion-=arc_dispersion(arc_box[right_index],ω)
            deleteat!(arc_box, right_index)
            insert!(arc_box, right_index, new_arc_r)
            # diagram.total_dispersion+=arc_dispersion(new_arc_r,ω)
        end

        if new_arc_l.index_out-new_arc_l.index_in == 2
            line_box[new_arc_l.index_in+1].covered=true
        end

        if new_arc_r.index_out-new_arc_r.index_in == 2
            line_box[new_arc_r.index_in+1].covered=true
        end
        sign=diagram.sign_box[line_index]

        shift_3!(diagram,sign[1],line_index-1,right_index,!right_open,m,μ,order)
        shift_3!(diagram,sign[2],line_index,left_index,!left_open,m,μ,order)
        # update_arcp!(diagram,order,right_index,!right_open,m,μ)
        # update_arcp!(diagram,order,left_index,!left_open,m,μ)#update_arcp_2!
        # shift!(diagram,sign[1],line_index-1,right_index,!right_open,m,μ)
        # shift!(diagram,sign[2],line_index,left_index,!left_open,m,μ)


        return true
    end
end

function insert_prob(arc::Arc,open_arc::Bool,ω::Int64,m::Int64,τ::Float64)
    q=arc.q
    period=arc.period
    if open_arc
        arc_T=τ*(1+period[2]-period[1])
    else
        arc_T=τ*(period[2]-period[1])
    end

    return ω*exp(-ω*(arc_T))*exp(-norm(q)^2/(2m)*arc_T)/(2pi*m/arc_T)^0.5#/norm(q)^2*1.0/(2pi)
end
