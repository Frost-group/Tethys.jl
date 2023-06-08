function insert_arc!(diagram::Diagram,order::Int64,m::Int64,μ::Float64,ω::Int64,α_squared::Float64)

    τ=diagram.τ
    cross_over=false

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

    if arc_T > τ
        return false
    end

    line.period[1]=τ_1/τ
    arc_T=τ_2-τ_1
    # q=MVector{3}(rand(Normal(0,sqrt(m/arc_T)),3))
    phi = rand(Uniform(0,pi*2))
    costheta = rand(Uniform(-1,1))
    theta = acos(costheta)
    x = sin(theta)*cos(phi)
    y = sin(theta)*sin(phi)
    z = cos(theta)
    q = abs(rand(Normal(0,sqrt(m/arc_T)))).*[x,y,z]

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
    else
        cross_over=true
        τ_2=τ_2-τ
        for i in index_in:2order+1
            line_tem=line_box[i]
            total_dis+=dispersion(line_tem, m, μ)*τ
            line_tem.k-=q
            line_tem.index+=1
            total_dis-=dispersion(line_tem, m, μ)*τ
        end

        for i in 1:index_in
            line_tem=line_box[i]

            if line_tem.period[2]<τ_2/τ
                total_dis+=dispersion(line_tem, m, μ)*τ
                line_tem.k-=q
                total_dis-=dispersion(line_tem, m, μ)*τ
                continue
            elseif i != index_in
                k_out=deepcopy(line_tem.k)
                τ_R_2=deepcopy(line_tem.period[2]*τ)
                line_tem.period[2]=τ_2/τ
                total_dis+=dispersion(line_tem, m, μ)*τ
                line_tem.k-=q
                total_dis-=dispersion(line_tem, m, μ)*τ
                index_out=i+1
                break
            else
                k_out=deepcopy(line_tem.k)
                line_tem.period[1]=τ_L/τ
                line_tem.period[2]=τ_2/τ
                total_dis-=dispersion(line_tem, m, μ)*τ
                line_tem.k+=q
                total_dis+=dispersion(line_tem, m, μ)*τ
                index_out=i+1
                line_tem.period[1]=τ_2/τ
                line_tem.period[2]=τ_1/τ
            end
        end

        index_in=index_in+1
    end


    new_arc=Arc(q,[τ_1,τ_2]/τ,ω,index_in,index_out)

    p_x_y=diagram.p_ins*ω/(1-exp(-ω*(τ)))/(2order+1)/(τ_R-τ_L)
    p_x_y*=1/(2pi)*exp(-norm(q)^2/(2m)*arc_T)/(2pi*m/arc_T)^0.5#*1.0/norm(q)^2
    p_y_x=diagram.p_rem/(order+1)
    r=α_squared*p_y_x/(exp(total_dis)*p_x_y*(2*pi)^3)#*norm(q)^2
    # coef_old=-diagram.dispersion/τ
    # coef_new=-(diagram.dispersion-total_dis-arc_T*ω)/τ
    # r*=1/(2order+2)/(2order+1)*(coef_new/coef_old)^(2order+1)*coef_new^2

    if r<rand()
        if !cross_over
            line_box[index_in].period[1]=τ_L/τ
            line_box[index_out-2].period[2]=τ_R_2/τ
            for i in index_in:index_out-2
                line_tem=line_box[i]
                line_tem.index=i
                line_tem.k+=q
            end
        else
            if index_in!=index_out
                line_box[index_in-1].period[1]=τ_L/τ
                line_box[index_out-1].period[2]=τ_R_2/τ
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
                line_box[index_in-1].period[1]=τ_L/τ
                line_box[index_in-1].period[2]=τ_R/τ
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
        diagram.order+=1
        diagram.dispersion-=total_dis
        diagram.dispersion-=arc_T*ω 

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
            line_tem=Line(k_in,[τ_L,τ_1]/τ, index_in, false)
            insert!(line_box, index_in, line_tem)
            line_tem=Line(k_out,[τ_2,τ_R_2]/τ, index_out, false)
            insert!(line_box, index_out, line_tem)
        else
            if index_out!=index_in
                line_tem=Line(k_out,[τ_2,τ_R_2]/τ, index_out, false)
                insert!(line_box, index_out, line_tem)
                line_tem=Line(k_in,[τ_L,τ_1]/τ, index_in, false)
                insert!(line_box, index_in, line_tem)
            else
                line_tem=Line(k_out,[τ_L,τ_2]/τ, index_in-1, false)
                insert!(line_box, index_out-1, line_tem)
                line_tem=Line(k_out,[τ_1,τ_R]/τ, index_in+1, false)
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

        # for i in 1:length(line_box)
        #     line_box[i].index=i
        #     total_dispersion+=dispersion(line_box[i], m, μ)
        # end

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

        # diagram.total_dispersion = total_dispersion
        return true
    end
end

function remove_arc!(diagram::Diagram,order::Int64,m::Int64,μ::Float64,ω::Int64,α_squared::Float64)

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
    line_in=line_box[index_in]
    line_out=line_box[index_out]


    τ_L=line_in.period[1]*τ

    if index_out-index_in==2
        τ_R=line_out.period[2]*τ
    elseif index_out==index_in
        τ_L=line_box[index_in-1].period[1]*τ
        τ_R=line_box[index_in+1].period[2]*τ
    else
        τ_R=line_box[index_in+1].period[2]*τ
    end

    τ_R_2=line_out.period[2]*τ
    
    τ_1=arc.period[1]*τ
    τ_2=arc.period[2]*τ
    arc_T=abs(offset_τ-abs(τ_2-τ_1))
    total_dis=0
    w_x=1
    w_y=exp(-ω*arc_T)*α_squared/(2*pi)^3/norm(q)^2

    open_arc_range = [collect(1:index_out-1); collect(index_in+1:2order+1)]

    if closed_arc
        for i in index_in+1:index_out-1
            line_tem=line_box[i]
            total_dis+=dispersion(line_tem, m, μ)*τ
            line_tem.index-=1
            line_tem.k+=q
            total_dis-=dispersion(line_tem, m, μ)*τ
        end
    else
        for i in open_arc_range
            line_tem=line_box[i]
            total_dis+=dispersion(line_tem, m, μ)*τ
            line_tem.k+=q
            total_dis-=dispersion(line_tem, m, μ)*τ
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

    p_x_y=diagram.p_ins*ω/(1-exp(-ω*(τ)))/(2order-1)/(τ_R-τ_L)
    p_x_y*=1.0/(2pi)*exp(-norm(q)^2/(2m)*arc_T)/(2pi*m/arc_T)^0.5#exp(-norm(q)^2/(2m)*arc_T)/(2pi*m/arc_T)^1.5

    p_y_x=diagram.p_rem/order
    r=((2*pi)^3*p_x_y)/(exp(total_dis)*p_y_x*α_squared)#*norm(q)^2
    # coef_old=-diagram.dispersion/τ
    # coef_new=-(diagram.dispersion-total_dis+arc_T*ω)/τ
    # r*=(2order)*(2order-1)*(coef_new/coef_old)^(2order-1)/coef_old^2

    if r<rand()
        if closed_arc
            for i in index_in+1:index_out-1
                line_tem=line_box[i]
                line_tem.index+=1
                line_tem.k-=q
            end
        else
            for i in open_arc_range
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
        sign_box=diagram.sign_box
        diagram.order-=1
        diagram.dispersion-=total_dis
        diagram.dispersion+=arc_T*ω 

        if closed_arc
            deleteat!(arc_box, index)
        else
            deleteat!(end_arc_box, index-arc_box_length)
        end

        if closed_arc
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

        else
            if index_out == index_in
                line_tem=line_box[index_in]
                line_tem.period[1]=τ_L/τ
                line_tem.period[2]=τ_R/τ
                line_tem.covered=false
                sign_to_add=[sign_box[index_out-1][1],sign_box[index_in+1][2]]
                deleteat!(sign_box, index_out-1:index_in+1)
                insert!(sign_box, index_out-1, sign_to_add)

            else
                line_tem=line_box[index_in]
                line_tem.period[1]=τ_L/τ
                line_tem.period[2]=τ_R/τ
                sign_to_add=[sign_box[index_in][1],sign_box[index_in+1][2]]
                deleteat!(sign_box,index_in:index_in+1)
                insert!(sign_box, index_in, sign_to_add)

                line_tem=line_box[index_out]
                line_tem.period[1]=line_box[index_out-1].period[1]
                sign_to_add=[sign_box[index_out-1][1],sign_box[index_out][2]]
                deleteat!(sign_box,index_out-1:index_out)
                insert!(sign_box, index_out-1, sign_to_add)

            end
            deleteat!(line_box, [index_out-1, index_in+1])
        end

        # total_dispersion = 0.0

        for i in 1:length(line_box)
            line_box[i].index=i
            # total_dispersion+=dispersion(line_box[i], m, μ)
        end

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

        # diagram.total_dispersion = total_dispersion
        return true
    end
end

function resample_arc!(diagram::Diagram,order::Int64,m::Int64,μ::Float64,ω::Int64,α_squared::Float64)

    # println("start")
    # println(diagram.dispersion)
    if order-1<0
        return false
    end
    
    arc_box=diagram.arc_box
    end_arc_box = diagram.end_arc_box
    arc_box_length = length(arc_box)
    index=rand(1:order)
    # println(index)
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
    index_out_r=arc.index_out
    # println(index_out_r)
    q=arc.q


    line_box=diagram.line_box
    line_in=line_box[index_in]
    line_out=line_box[index_out_r]

    τ_1=arc.period[1]*τ
    τ_2=arc.period[2]*τ
    arc_T=abs(offset_τ-abs(τ_2-τ_1))

    τ_2_in=τ_1-log(rand())/ω
    arc_T_in=τ_2_in-τ_1
    if arc_T_in > τ
        println("overlap")
        return false
    end

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


    cross_over=false
    phi = rand(Uniform(0,pi*2))
    costheta = rand(Uniform(-1,1))
    theta = acos(costheta)
    x = sin(theta)*cos(phi)
    y = sin(theta)*sin(phi)
    z = cos(theta)
    q_in = abs(rand(Normal(0,sqrt(m/arc_T_in)))).*[x,y,z]

    τ_R_2=0.0
    index_out_in=0
    k_out=0

    if τ_2_in<τ
        #not set covered yet
        for i in index_in+1:2order+1
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
                line_tem.index+=1
                total_dis+=dispersion(line_tem, m, μ)*τ
                index_out_in=i+1
                break
            end
        end
    else
        cross_over=true
        τ_2_in=τ_2_in-τ
        for i in index_in+1:2order+1
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
        # index_in=index_in+1
    end


    new_arc=Arc(q_in,[τ_1,τ_2_in]/τ,ω,index_in,index_out_in)

    p_x_y=exp(-norm(q)^2/(2m)*arc_T)/(2pi*m/arc_T)^0.5#*exp(-ω*arc_T)
    p_x_y/=exp(-norm(q_in)^2/(2m)*arc_T_in)/(2pi*m/arc_T_in)^0.5#*exp(-ω*arc_T_in)
    r=exp(total_dis)*p_x_y
    # println(τ_2_in/τ)
    # println(index_out_in)
    # println(index_out_r)
    if r<rand()

        if !cross_over
            line_box[index_out_in-1].period[2]=τ_R_2
            for i in index_in+1:index_out_in-1
                line_tem=line_box[i]
                line_tem.index=i
                line_tem.k+=q_in
            end
        else
            line_box[index_out_in-1].period[2]=τ_R_2
            for i in index_in+1:2order+1
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
        # println("closed_arc")
        # println(closed_arc)
        # println("cross_over")
        # println(cross_over)
        sign_box=diagram.sign_box
        diagram.dispersion+=total_dis
        diagram.dispersion+=(arc_T-arc_T_in)*ω 
        # arc_box=diagram.arc_box
        # end_arc_box = diagram.end_arc_box
        if closed_arc
            deleteat!(arc_box, index)
        else
            deleteat!(end_arc_box, index-arc_box_length)
        end
        # println(diagram.arc_box)
        # println(q_in)

        if index_out_r == index_out_in
            # println("s1")
            line_tem=line_box[index_out_r]
            line_tem.period[1]=τ_2_in/τ
        elseif index_out_r == index_out_in-1
            # println("s2")
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
            if closed_arc
                if index_out_in>index_out_r
                    # println("s3")
                    index_out_in-=1
                    new_arc.index_out=index_out_in
                    insert!(line_box, index_out_in, line_tem_3)
                    sign_to_add=[[sign_box[index_out_in-1][1],1],[1,sign_box[index_out_in-1][2]]]
                    deleteat!(sign_box,index_out_in-1)

                    for i in 1:2
                        insert!(sign_box, index_out_in-1, sign_to_add[3-i])
                    end
                elseif index_out_in<index_out_r && !cross_over
                    # println("s4")
                    insert!(line_box, index_out_in, line_tem_3)
                    sign_to_add=[[sign_box[index_out_in-1][1],1],[1,sign_box[index_out_in-1][2]]]
                    deleteat!(sign_box,index_out_in-1)

                    for i in 1:2
                        insert!(sign_box, index_out_in-1, sign_to_add[3-i])
                    end
                else
                    # println("s5")
                    new_arc.index_in+=1
                    insert!(line_box, index_out_in, line_tem_3)
                    sign_to_add=[[sign_box[index_out_in-1][1],1],[1,sign_box[index_out_in-1][2]]]
                    deleteat!(sign_box,index_out_in-1)

                    for i in 1:2
                        insert!(sign_box, index_out_in-1, sign_to_add[3-i])
                    end
                end
            else
                if index_out_in>index_out_r && cross_over
                    index_out_in-=1
                    new_arc.index_out=index_out_in
                    insert!(line_box, index_out_in, line_tem_3)
                    sign_to_add=[[sign_box[index_out_in-1][1],1],[1,sign_box[index_out_in-1][2]]]
                    deleteat!(sign_box,index_out_in-1)

                    for i in 1:2
                        insert!(sign_box, index_out_in-1, sign_to_add[3-i])
                    end
                elseif index_out_in>index_out_r && !cross_over
                    index_out_in-=1
                    new_arc.index_out=index_out_in
                    new_arc.index_in-=1
                    insert!(line_box, index_out_in, line_tem_3)
                    sign_to_add=[[sign_box[index_out_in-1][1],1],[1,sign_box[index_out_in-1][2]]]
                    deleteat!(sign_box,index_out_in-1)

                    for i in 1:2
                        insert!(sign_box, index_out_in-1, sign_to_add[3-i])
                    end
                else
                    # new_arc.index_in+=1
                    insert!(line_box, index_out_in, line_tem_3)
                    sign_to_add=[[sign_box[index_out_in-1][1],1],[1,sign_box[index_out_in-1][2]]]
                    deleteat!(sign_box,index_out_in-1)

                    for i in 1:2
                        insert!(sign_box, index_out_in-1, sign_to_add[3-i])
                    end
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
        # println("insert")
        if !cross_over
            push!(diagram.arc_box,new_arc)
        else
            push!(diagram.end_arc_box,new_arc)
            diagram.component+=1
        end

        return true
    end
end


# dia=diagram.dispersion
# new_dia= dia-total_dis+arc_T*ω
# if dia !=0 && new_dia !=0
#     factor=(2order)*(2order-1)*(new_dia/dia)^(2order-1)
#     factor/=(dia/τ)^2
# elseif dia == 0 && new_dia !=0
#     factor=(-new_dia/τ)^(2order-1)/factorial(2order-2)
# else
#     factor=factorial(2order)/(-dia/τ)^(2order+1)
# end

# dia=diagram.dispersion
# new_dia= dia-total_dis-arc_T*ω
# if dia !=0 && new_dia !=0
#     factor=(2order+2)*(2order+1)*(dia/new_dia)^(2order+1)
#     factor/=(new_dia/τ)^2
# elseif dia == 0 && new_dia !=0
#     factor=factorial(2order+2)/(-new_dia/τ)^(2order+3)
# else
#     factor=(-dia/τ)^(2order+1)/factorial(2order)
# end

using SpecialFunctions

2pi^2m/(norm(mean_k)*arc_T)*erfi(sqrt(norm(mean_k)^2*arc_T/(2m)))

2pi*sqrt(2pim/arc_T)