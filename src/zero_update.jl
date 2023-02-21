include("Diagram.jl")
using Base.Threads
using StructArrays


function zero_insert_arc!(diagram::Diagram,order::Int64,m::Int64,μ::Float64,ω::Int64,α_squared::Float64)

    τ=diagram.τ
    component=diagram.component

    if order+1>diagram.max_order
        return false
    end

    line_box=diagram.line_box
    index_in=rand(1:length(line_box))
    line=line_box[index_in]
    τ_L=deepcopy(line.period[1])*τ
    τ_R=deepcopy(line.period[2])*τ
    k_in=deepcopy(line.k)

    τ_1=rand(Uniform(τ_L,τ_R))
    τ_2=τ_1-log(rand())/ω

    arc_T=τ_2-τ_1

    if τ_2 > τ
        return false
    end

    line.period[1]=τ_1/τ
    arc_T=τ_2-τ_1
    q=MVector{3}(rand(Normal(0,sqrt(m/arc_T)),3))
    
    w_x=1.0
    w_y=1.0
    total_dis=0.0
    τ_R_2=0.0
    index_out=0
    k_out=0

    #not set covered yet
    for i in index_in:2order+1
        line_tem=line_box[i]
        if line_tem.period[2]<τ_2/τ
            total_dis+=dispersion(line_tem, m, μ)*τ
            line_tem.k-=q
            line_tem.index+=1
            total_dis-=dispersion(line_tem, m, μ)*τ
            continue
        else
            k_out=deepcopy(line_tem.k)
            τ_R_2=deepcopy(line_tem.period[2]*τ)
            line_tem.period[2]=τ_2/τ
            total_dis+=dispersion(line_tem, m, μ)*τ
            line_tem.k-=q
            line_tem.index+=1
            total_dis-=dispersion(line_tem, m, μ)*τ
            index_out=i+2
            break
        end
    end


    new_arc=Arc(q,[τ_1,τ_2]/τ,ω,index_in,index_out)

    p_x_y=diagram.p_ins/(2order+1)/(τ_R-τ_L)
    p_x_y*=exp(-norm(q)^2/(2m)*arc_T)/(2pi*m/arc_T)^1.5
    p_y_x=diagram.p_rem/(order+1-component)
    r=α_squared*p_y_x/(exp(total_dis)*p_x_y*(2*pi)^3*norm(q)^2)
    # coef_old=-diagram.dispersion/τ
    # coef_new=-(diagram.dispersion-total_dis-arc_T*ω)/τ
    # r*=1/(2order+2)/(2order+1)*(coef_new/coef_old)^(2order+1)*coef_new^2

    if r<rand()

        line_box[index_in].period[1]=τ_L/τ
        line_box[index_out-2].period[2]=τ_R_2/τ
        for i in index_in:index_out-2
            line_tem=line_box[i]
            line_tem.index=i
            line_tem.k+=q
        end

        return false
    else
        diagram.order+=1
        diagram.dispersion-=total_dis
        diagram.dispersion-=arc_T*ω 

        if index_out-index_in==2
            line_box[index_in].covered=true
        else
            line_box[index_in].covered=false
            line_box[index_out-2].covered=false
        end

        line_tem=Line(k_in,[τ_L,τ_1]/τ, index_in, false)
        insert!(line_box, index_in, line_tem)
        line_tem=Line(k_out,[τ_2,τ_R_2]/τ, index_out, false)
        insert!(line_box, index_out, line_tem)

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
        end

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

        push!(diagram.arc_box,new_arc)

        # diagram.total_dispersion = total_dispersion
        return true
    end
end

function zero_remove_arc!(diagram::Diagram,order::Int64,m::Int64,μ::Float64,ω::Int64,α_squared::Float64)

    if order-1<0
        return false
    end
    
    arc_box=diagram.arc_box
    component=diagram.component
    arc_box_length = length(arc_box)
    index=rand(1:order-component)
    τ=diagram.τ
    arc=arc_box[index]
        
    index_in=arc.index_in
    index_out=arc.index_out
    q=arc.q

    line_box=diagram.line_box
    line_in=line_box[index_in]
    line_out=line_box[index_out]


    τ_L=line_in.period[1]*τ

    if index_out-index_in==2
        τ_R=line_out.period[2]*τ
    else
        τ_R=line_box[index_in+1].period[2]*τ
    end

    τ_R_2=line_out.period[2]*τ
    
    τ_1=arc.period[1]*τ
    τ_2=arc.period[2]*τ
    arc_T=abs(τ_2-τ_1)
    total_dis=0
    w_x=1
    w_y=exp(-ω*arc_T)*α_squared/(2*pi)^3/norm(q)^2

    for i in index_in+1:index_out-1
        line_tem=line_box[i]
        total_dis+=dispersion(line_tem, m, μ)*τ
        line_tem.index-=1
        line_tem.k+=q
        total_dis-=dispersion(line_tem, m, μ)*τ
    end

    p_x_y=diagram.p_ins/(2order-1)/(τ_R-τ_L)
    p_x_y*=exp(-norm(q)^2/(2m)*arc_T)/(2pi*m/arc_T)^1.5

    p_y_x=diagram.p_rem/(order-component)
    r=((2*pi)^3*p_x_y*norm(q)^2)/(exp(total_dis)*p_y_x*α_squared)
    # coef_old=-diagram.dispersion/τ
    # coef_new=-(diagram.dispersion-total_dis+arc_T*ω)/τ
    # r*=(2order)*(2order-1)*(coef_new/coef_old)^(2order-1)/coef_old^2

    if r<rand()
        for i in index_in+1:index_out-1
            line_tem=line_box[i]
            line_tem.index+=1
            line_tem.k-=q
        end
        return false
    else
        sign_box=diagram.sign_box
        diagram.order-=1
        diagram.dispersion-=total_dis
        diagram.dispersion+=arc_T*ω 


        deleteat!(arc_box, index)

        if index_out-index_in==2
            line_tem=line_box[index_in+1]
            line_tem.period[1]=τ_L/τ
            line_tem.period[2]=τ_R/τ
            line_tem.covered=false
            sign_to_add=[sign_box[index_in][1],sign_box[index_out][2]]
            deleteat!(sign_box, index_in:index_out)
            insert!(sign_box, index_in, sign_to_add)

        else
            line_tem=line_box[index_in+1]
            line_tem.period[1]=τ_L/τ
            sign_to_add=[sign_box[index_in][1],sign_box[index_in+1][2]]
            deleteat!(sign_box, index_in:index_in+1)
            insert!(sign_box, index_in, sign_to_add)


            line_tem=line_box[index_out-1]
            line_tem.period[2]=τ_R_2/τ
            sign_to_add=[sign_box[index_out-2][1],sign_box[index_out-1][2]]
            deleteat!(sign_box, index_out-2:index_out-1)
            insert!(sign_box, index_out-2, sign_to_add)

        end

        deleteat!(line_box, [index_in,index_out])

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

        return true
    end
end


function zero_swap_arc!(diagram::Diagram)

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
    arc_box_length = length(arc_box)
    left_check=false
    right_check=false
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

    new_line=Line(chosen_line.k-sign[1]*q1+sign[2]*q2 ,chosen_line.period, line_index,false)
    w_x=green_zero(chosen_line, m, μ)*phonon_propagator(arc_l)*phonon_propagator(arc_r)
    w_y=green_zero(new_line, m, μ)*phonon_propagator(new_arc_l)*phonon_propagator(new_arc_r)

    r=(w_y/w_x)^diagram.τ

    if r<rand()
        return false
    else
        diagram.dispersion+=log(r)
        deleteat!(line_box, line_index)
        insert!(line_box, line_index, new_line)
        # diagram.total_dispersion+=dispersion(new_line, m, μ)
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
