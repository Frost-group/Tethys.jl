include("Diagram.jl")
using Base.Threads
using StructArrays
using SpecialFunctions
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

function insert_arc!(diagram::Diagram,order::Int64,m::Int64,μ::Float64,ω::Int64,α_squared::Float64)

    τ=diagram.τ
    cross_over=false

    if order+1>diagram.max_order
        return false
    end

    line_box=diagram.line_box
    index_in=rand(1:length(line_box))
    line=line_box[index_in]
    # τ_1=rand(Uniform(0,1))
    # i=1
    # while τ_1>line_box[i].period[2]
    #     i+=1
    # end
    # index_in=i
    # line=line_box[index_in]

    τ_L=deepcopy(line.period[1])*τ
    τ_R=deepcopy(line.period[2])*τ
    k_in=deepcopy(line.k)
    # τ_1*=τ
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
    # costheta = rand(Uniform(-1,1))
    # theta = acos(costheta)
    # x = sin(theta)*cos(phi)
    # y = sin(theta)*sin(phi)
    # z = cos(theta)
    q = abs(rand(Normal(0,sqrt(m/arc_T))))#.*[x,y,z]

    w_x=1.0
    w_y=1.0
    total_dis=0.0
    τ_R_2=0.0
    index_out=0
    k_out=0
    mean_k=0

    if τ_2<τ
        #not set covered yet
        for i in index_in:2order+1
            line_tem=line_box[i]
            if line_tem.period[2]<τ_2/τ
                k=line_tem.k
                period=line_tem.period
                mean_k=mean_k.+k*(period[2]-period[1])
                continue
            else
                k=line_tem.k
                period=line_tem.period
                k_out=deepcopy(line_tem.k)
                τ_R_2=deepcopy(line_tem.period[2]*τ)
                line_tem.period[2]=τ_2/τ
                mean_k=mean_k.+k*(period[2]-period[1])
                break
            end
        end
    else
        cross_over=true
        τ_2=τ_2-τ
        for i in index_in:2order+1
            line_tem=line_box[i]
            k=line_tem.k
            period=line_tem.period
            mean_k=mean_k.+k*(period[2]-period[1])
        end
        for i in 1:index_in
            line_tem=line_box[i]
            if line_tem.period[2]<τ_2/τ
                k=line_tem.k
                period=line_tem.period
                mean_k=mean_k.+k*(period[2]-period[1])
                continue
            elseif i != index_in
                k_out=deepcopy(line_tem.k)
                τ_R_2=deepcopy(line_tem.period[2]*τ)
                line_tem.period[2]=τ_2/τ
                k=line_tem.k
                period=line_tem.period
                mean_k=mean_k.+k*(period[2]-period[1])
                break
            else
                k=line_tem.k
                period=line_tem.period
                mean_k=mean_k.+k*(τ_2/τ-τ_L/τ)
            end
        end
    end

    
    if norm(mean_k)<=1e-10
        costheta = rand(Uniform(-1,1))
        theta = acos(costheta)
        prob=1.0/2.
    else
        coe=q*norm(mean_k)*τ/m
        costheta = 1+log(1-rand()*(1-exp(-2coe)))/coe
        # println(norm(mean_k))
        theta = acos(costheta)
        prob=coe/(exp(coe*(1-costheta))-exp(coe*(-1-costheta)))
        # println(sinh(coe)/coe)
    end

    x = sin(theta)*cos(phi)
    y = sin(theta)*sin(phi)
    z = cos(theta)
    direct=[x,y,z]

    if norm(mean_k) !=0 
        direct=rotate(mean_k)*direct
        # println(mean_k/norm(mean_k))
        # println(rotate(mean_k)*[0,0,1])
    end
    
    # println(dot(mean_k,direct)/norm(mean_k))
    # println(cos(theta))
    q=q*direct

    if !cross_over
        #not set covered yet
        for i in index_in:2order+1
            line_tem=line_box[i]
            if line_tem.period[2]<τ_2/τ
                total_dis+=dispersion(line_tem, m, μ)*τ
                line_tem.k-=q
                total_dis-=dispersion(line_tem, m, μ)*τ
                continue
            else
                total_dis+=dispersion(line_tem, m, μ)*τ
                line_tem.k-=q
                total_dis-=dispersion(line_tem, m, μ)*τ
                index_out=i+2
                break
            end
        end
    else
        for i in index_in:2order+1
            line_tem=line_box[i]
            total_dis+=dispersion(line_tem, m, μ)*τ
            line_tem.k-=q
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
    p_x_y*=2/(2pi)*exp(-norm(q)^2/(2m)*arc_T)/(2pi*m/arc_T)^0.5*prob#*1.0/norm(q)^2
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

function rotate(k)
    e=k/norm(k)
    x=acos(e[3])
    u=[-e[2]/sin(x),e[1]/sin(x),0]
    array_1=[cos(x)+u[1]^2*(1-cos(x)), u[1]*u[2]*(1-cos(x)), -u[2]*sin(x)]
    array_2=[u[1]*u[2]*(1-cos(x)), cos(x)+u[2]^2*(1-cos(x)), u[1]*sin(x)]
    array_3=[u[2]*sin(x), -u[1]*sin(x), cos(x)]

    if x==0
        return[[1,0,0] [0,1,0] [0,0,1]]
    end
    return [ array_1 array_2 array_3]
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
    mean_k=0
    
    if closed_arc
        for i in index_in+1:index_out-1
            line_tem=line_box[i]
            total_dis+=dispersion(line_tem, m, μ)*τ
            line_tem.index-=1
            line_tem.k+=q
            k=line_tem.k
            period=line_tem.period
            mean_k=mean_k.+k*(period[2]-period[1])
            total_dis-=dispersion(line_tem, m, μ)*τ
        end
    else
        for i in open_arc_range
            line_tem=line_box[i]
            total_dis+=dispersion(line_tem, m, μ)*τ
            line_tem.k+=q
            k=line_tem.k
            period=line_tem.period
            mean_k=mean_k.+k*(period[2]-period[1])
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

    if norm(mean_k)<=1e-10
        prob=1.0/2.
    else
        costheta=dot(mean_k,q)/(norm(mean_k)*norm(q))
        coe=norm(q)*norm(mean_k)*τ/m
        prob=coe/(exp(coe*(1-costheta))-exp(coe*(-1-costheta)))
    end

    p_x_y=diagram.p_ins*ω/(1-exp(-ω*(τ)))/(2order-1)/(τ_R-τ_L)
    p_x_y*=2.0/(2pi)*exp(-norm(q)^2/(2m)*arc_T)/(2pi*m/arc_T)^0.5*prob#exp(-norm(q)^2/(2m)*arc_T)/(2pi*m/arc_T)^1.5

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

    if τ_new<9.0
        return false
    end

    line_end.period[2]=τ_new
    line_box[2*order+1]=line_end
    diagram.τ=τ_new

    return true

end


function scale!(diagram::Diagram, order::Int64,m::Int64,μ::Float64,ω::Int64, samples::Int64)

    line_box=diagram.line_box
    arc_box=diagram.arc_box
    end_arc_box=diagram.end_arc_box
    #line_box_length = length(line_box)
    #arc_box_length = length(arc_box)
    #end_arc_box_length = length(end_arc_box)
    τ=diagram.τ
    record_τ=diagram.record_τ
    total_dis=diagram.dispersion
    
    coef=-τ/total_dis
    #println(coef)
    τ_new=rand(Gamma(2*order+1,coef),samples)

    if τ_new[1]>diagram.max_τ
        return false
    end

    diagram.record_τ=τ_new

    # if τ_new[1]<5.0
    #     return false
    # end
    
    diagram.dispersion*=τ_new[1]/τ
    diagram.τ=τ_new[1]

    return true

end

function set_μ!(diagram::Diagram,new_μ::Float64)
    old_μ=copy(diagram.μ)
    total_dis=diagram.dispersion
    τ=diagram.τ

    diagram.dispersion+=old_μ*τ
    diagram.dispersion-=new_μ*τ
    diagram.μ=new_μ

    return diagram
end

function set_τ!(diagram::Diagram,new_τ::Float64)

    τ=diagram.τ
    diagram.dispersion/=τ
    diagram.dispersion*=new_τ
    diagram.τ=new_τ

    return diagram
end

function energy(diagram::Diagram)
    total_dis=-diagram.dispersion
    τ=diagram.τ
    order=diagram.order
    #return (total_dis-2*order+2*length(diagram.end_arc_box))/τ
    return (total_dis-2*order)/τ
end

function mass_estimator(diagram::Diagram)
    line_box=diagram.line_box
    p_squared = [0.0, 0.0, 0.0]
    τ=diagram.τ
    m=diagram.mass
    μ=diagram.μ
    for line in line_box
        p_squared .+= line.k*(line.period[2]-period[1])#p_dispersion(line, m, μ) .*2
    end
    #return 1 / (1 - (norm(p_squared)^2)*τ/(3))
    return (norm(p_squared)^2)
end

function total_dis_check(diagram::Diagram)

    m=diagram.mass
    ω=diagram.ω
    μ=diagram.μ
    τ=diagram.τ

    total_dis=0

    for line in diagram.line_box
        total_dis+=dispersion(line, m, μ)*τ
    end

    for arc in diagram.arc_box
        total_dis+=arc_dispersion(arc, ω)*τ
    end

    for arc in diagram.end_arc_box
        total_dis+=-ω*τ+arc_dispersion(arc, ω)*τ
    end

    return total_dis
end

function resample_arc!(diagram::Diagram,order::Int64,m::Int64,μ::Float64,ω::Int64,α_squared::Float64)

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

    open_arc_range = [collect(1:index_out_r-1); collect(index_in+1:2order+1)]
    mean_k_r=0

    if closed_arc
        for i in index_in+1:index_out_r-1
            line_tem=line_box[i]
            total_dis-=dispersion(line_tem, m, μ)*τ
            line_tem.index-=1
            line_tem.k+=q
            k=line_tem.k
            period=line_tem.period
            mean_k_r=mean_k_r.+k*(period[2]-period[1])
            total_dis+=dispersion(line_tem, m, μ)*τ
        end
    else
        for i in open_arc_range
            line_tem=line_box[i]
            total_dis-=dispersion(line_tem, m, μ)*τ
            line_tem.k+=q
            k=line_tem.k
            period=line_tem.period
            mean_k_r=mean_k_r.+k*(period[2]-period[1])
            total_dis+=dispersion(line_tem, m, μ)*τ
        end        
    end

    if norm(mean_k_r)<=1e-10
        prob_r=1.0/2.
    else
        costheta=dot(mean_k_r,q)/(norm(mean_k_r)*norm(q))
        coe=norm(q)*norm(mean_k_r)*τ/m
        prob_r=coe/(exp(coe*(1-costheta))-exp(coe*(-1-costheta)))
    end

    #insert
    cross_over=false
    phi = rand(Uniform(0,pi*2))
    q_in = abs(rand(Normal(0,sqrt(m/arc_T_in))))

    τ_R_2=0.0
    index_out_in=0
    k_out=0

    mean_k_i=0

    if τ_2_in<τ
        #not set covered yet
        for i in index_in+1:2order+1
            line_tem=line_box[i]
            if line_tem.period[2]<τ_2_in/τ
                k=line_tem.k
                period=line_tem.period
                mean_k_i=mean_k_i.+k*(period[2]-period[1])
                continue
            else
                k=line_tem.k
                period=line_tem.period
                k_out=deepcopy(line_tem.k)
                τ_R_2=deepcopy(line_tem.period[2])
                line_tem.period[2]=τ_2_in/τ
                mean_k_i=mean_k_i.+k*(period[2]-period[1])
                break
            end
        end
    else
        cross_over=true
        τ_2_in=τ_2_in-τ
        for i in index_in+1:2order+1
            line_tem=line_box[i]
            k=line_tem.k
            period=line_tem.period
            mean_k_i=mean_k_i.+k*(period[2]-period[1])
        end
        for i in 1:index_in
            line_tem=line_box[i]
            if line_tem.period[2]<τ_2_in/τ
                k=line_tem.k
                period=line_tem.period
                mean_k_i=mean_k_i.+k*(period[2]-period[1])
                continue
            else
                k_out=deepcopy(line_tem.k)
                τ_R_2=deepcopy(line_tem.period[2])
                line_tem.period[2]=τ_2_in/τ
                k=line_tem.k
                period=line_tem.period
                mean_k_i=mean_k_i.+k*(period[2]-period[1])
                break
            end
        end
    end


    if norm(mean_k_i)<=1e-10
        costheta = rand(Uniform(-1,1))
        theta = acos(costheta)
        prob_i=1.0/2.
    else
        coe=q_in*norm(mean_k_i)*τ/m
        costheta = 1+log(1-rand()*(1-exp(-2coe)))/coe
        # println(norm(mean_k))
        theta = acos(costheta)
        prob_i=coe/(exp(coe*(1-costheta))-exp(coe*(-1-costheta)))
    end

    x = sin(theta)*cos(phi)
    y = sin(theta)*sin(phi)
    z = cos(theta)
    direct=[x,y,z]

    if norm(mean_k_i) !=0 
        direct=rotate(mean_k_i)*direct
    end

    q_in=q_in*direct

    if !cross_over
        #not set covered yet
        for i in index_in+1:2order+1
            line_tem=line_box[i]
            if line_tem.period[2]<τ_2_in/τ
                total_dis-=dispersion(line_tem, m, μ)*τ
                line_tem.k-=q_in
                total_dis+=dispersion(line_tem, m, μ)*τ
                continue
            else
                total_dis-=dispersion(line_tem, m, μ)*τ
                line_tem.k-=q_in
                total_dis+=dispersion(line_tem, m, μ)*τ
                index_out_in=i+1
                break
            end
        end
    else
        for i in index_in+1:2order+1
            line_tem=line_box[i]
            total_dis-=dispersion(line_tem, m, μ)*τ
            line_tem.k-=q_in
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
                total_dis-=dispersion(line_tem, m, μ)*τ
                line_tem.k-=q_in
                total_dis+=dispersion(line_tem, m, μ)*τ
                index_out_in=i+1
                break
            end
        end
    end

    new_arc=Arc(q_in,[τ_1,τ_2_in]/τ,ω,index_in,index_out_in)

    p_x_y=exp(-norm(q)^2/(2m)*arc_T)/(2pi*m/arc_T)^0.5#*exp(-ω*arc_T)
    p_x_y/=exp(-norm(q_in)^2/(2m)*arc_T_in)/(2pi*m/arc_T_in)^0.5#*exp(-ω*arc_T_in)
    r=exp(total_dis)*p_x_y*prob_r/prob_i
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

    if (chosen_line.period[2]-chosen_line.period[1])*ω*τ>-log(0.95)
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

        # shift_3!(diagram,sign[1],line_index-1,right_index,!right_open,m,μ,order)
        # shift_3!(diagram,sign[2],line_index,left_index,!left_open,m,μ,order)
        # update_arcp!(diagram,order,right_index,!right_open,m,μ)
        # update_arcp!(diagram,order,left_index,!left_open,m,μ)#update_arcp_2!
        # shift!(diagram,sign[1],line_index-1,right_index,!right_open,m,μ)
        # shift!(diagram,sign[2],line_index,left_index,!left_open,m,μ)


        return true
    end
end