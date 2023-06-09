include("Diagram.jl")
using Base.Threads
using StructArrays
using SpecialFunctions

function τ_update!(diagram::Diagram)

    """
    τ_update!(diagram::Diagram)

    This function is not in use.
    Extend the length of the 0th order Diagram object according to exponential distribution.
    """

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

    """
    p_update!(diagram::Diagram)

    This function is not in use.
    Sample the momentum along x-axis of the 0th order Diagram object according to Gaussian distribution.
    """

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

    """
    insert_arc!(diagram::Diagram,order::Int64,m::Int64,μ::Float64,ω::Int64,α_squared::Float64)

    Attempting to insert an Arc object into the Diagram object.

    # Arguments
    - `diagram::Diagram`: is the Diagram object to be modified.
    - `order::Int64`: is the number of Arc objects presented in the Diagram.
    - `m::Int64`: is the mass of the electron(assumed to be 1).
    - `μ::Float64`: is the energy offset for the Green's function(not in use).
    - `ω::Int64`: is the phonon frequency(assumed to be 1).
    - `α_squared::Float64`: is 2pi*sqrt(2) of coupling strength.
    """

    τ=diagram.τ
    cross_over=false # defines whether if the inserted Arc object is crossing the end of the Diagram

    # order checking
    if order+1>diagram.max_order
        return false
    end

    # uniform sampling for the line_in
    line_box=diagram.line_box
    index_in=rand(1:length(line_box)) # samples the line index attached to the left of the left vertex of the Arc
    line=line_box[index_in]

    τ_L=deepcopy(line.period[1])*τ
    τ_R=deepcopy(line.period[2])*τ
    k_in=deepcopy(line.k)
    τ_1=rand(Uniform(τ_L,τ_R)) # left vertex sampling (uniform distribution)
    τ_2=τ_1-log(rand())/ω # right vertex sampling (exponential distribution)

    arc_T=τ_2-τ_1

    # Arc length check
    if arc_T > τ
        return false
    end

    line.period[1]=τ_1/τ
    arc_T=τ_2-τ_1
    phi = rand(Uniform(0,pi*2)) # azimuth angle sampling for phonon momentum
    q = abs(rand(Normal(0,sqrt(m/arc_T)))) # phonon momentum magnitude sampling

    total_dis=0.0
    τ_R_2=0.0 # right vertex for the line_out attached to the Arc
    index_out=0 # line attached to the right of the right vertex
    k_out=0 # momentum for the line_out
    mean_k=0 # mean momentum of the Lines covered by the Arc

    # find out the region covered by the Arc
    if τ_2<τ
        # for the case not cross_over
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
        # for the case cross_over
        cross_over=true
        τ_2=τ_2-τ # modify the right vertex time

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

    # altitude angle sampling
    if norm(mean_k)<=1e-10
        # condition set for too small mean_k to avoid the error
        costheta = rand(Uniform(-1,1))
        theta = acos(costheta)
        prob=1.0/2.
    else
        coe=q*norm(mean_k)*τ/m
        costheta = 1+log(1-rand()*(1-exp(-2coe)))/coe
        theta = acos(costheta)
        prob=coe/(exp(coe*(1-costheta))-exp(coe*(-1-costheta)))
    end

    # transform  direction vector from spherical to Cartesian
    x = sin(theta)*cos(phi)
    y = sin(theta)*sin(phi)
    z = cos(theta)
    direct=[x,y,z]

    if norm(mean_k) !=0 
        direct=rotate(mean_k)*direct # rotate to correct direction
    end

    q=q*direct # final momentum

    # start to record the index_in & index out for Arc(after Arc is added) and the dispersion change
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


    new_arc=Arc(q,[τ_1,τ_2]/τ,ω,index_in,index_out) # Arc object to be inserted

    p_x_y=diagram.p_ins*ω/(1-exp(-ω*(τ)))/(2order+1)/(τ_R-τ_L)
    p_x_y*=2/(2pi)*exp(-norm(q)^2/(2m)*arc_T)/(2pi*m/arc_T)^0.5*prob
    p_y_x=diagram.p_rem/(order+1)
    r=α_squared*p_y_x/(exp(total_dis)*p_x_y*(2*pi)^3) # acceptance ratio

    if r<rand()
        # return back to original diagram
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
        # Arc insertation is accepted. Update information in diagram

        diagram.order+=1
        diagram.dispersion-=total_dis
        diagram.dispersion-=arc_T*ω 

        # check whether if the Arc covers a single line or not (not breaking by another Arc)
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

        # update the line_box in diagram (new Lines due to Arc insert)
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

        # diagram sign_box update
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

        # index_in & index_out update for Arcs in arc_box
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

        # add the Arc to the arc_box
        if !cross_over
            push!(diagram.arc_box,new_arc)
        else
            push!(diagram.end_arc_box,new_arc)
            diagram.component+=1 # record the crossing Arc
        end
        return true
    end
end

function rotate(k)
    """
    rotate(k)

    Return the rotation matrix to rotate z-axis to the mean-k direction

    # Arguments
    - `k::MVector{3,Float64}`: is mean_k for z-axis to map.
    """
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

    """
    remove_arc!(diagram::Diagram,order::Int64,m::Int64,μ::Float64,ω::Int64,α_squared::Float64)

    Attempting to remove an Arc object from the Diagram object.

    # Arguments
    - `diagram::Diagram`: is the Diagram object to be modified.
    - `order::Int64`: is the number of Arc objects presented in the Diagram.
    - `m::Int64`: is the mass of the electron(assumed to be 1).
    - `μ::Float64`: is the energy offset for the Green's function(not in use).
    - `ω::Int64`: is the phonon frequency(assumed to be 1).
    - `α_squared::Float64`: is 2pi*sqrt(2) of coupling strength.
    """

    # order checking
    if order-1<0
        return false
    end
    
    arc_box=diagram.arc_box
    end_arc_box = diagram.end_arc_box
    arc_box_length = length(arc_box)
    index=rand(1:order)
    offset_τ=0
    τ=diagram.τ

    # randomly select the Arc to be removed
    if index <= arc_box_length
        # for not crossing case
        arc=arc_box[index]
        closed_arc = true 
    else
        # for crossing case
        arc=end_arc_box[index-arc_box_length]
        closed_arc = false
        offset_τ=diagram.τ
    end

    # extract the Arc parameters
    index_in=arc.index_in
    index_out=arc.index_out
    q=arc.q

    line_box=diagram.line_box
    line_in=line_box[index_in]
    line_out=line_box[index_out]

    τ_L=line_in.period[1]*τ

    # work out the left vertex sampling range
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

    # to record the mean_k for phonon momentum direction sampling
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
    p_x_y*=2.0/(2pi)*exp(-norm(q)^2/(2m)*arc_T)/(2pi*m/arc_T)^0.5*prob
    p_y_x=diagram.p_rem/order
    r=((2*pi)^3*p_x_y)/(exp(total_dis)*p_y_x*α_squared) # acceptance ratio

    if r<rand()
    # return back to original diagram
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
        # Arc removal is accepted. Update information in diagram

        sign_box=diagram.sign_box
        diagram.order-=1
        diagram.dispersion-=total_dis
        diagram.dispersion+=arc_T*ω 

        # Arc delete
        if closed_arc
            deleteat!(arc_box, index)
        else
            deleteat!(end_arc_box, index-arc_box_length)
        end

        # update the covered property and the diagram sign box
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
            # Line object delete
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
            # Line object delete
            deleteat!(line_box, [index_out-1, index_in+1])
        end

        # set the line index to correct sequence (for debug)
        for i in 1:length(line_box)
            line_box[i].index=i
        end

        # index_in & index_out update for Arcs in arc_box
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

            diagram.component-=1 # record the crossing Arc
        end
        return true
    end
end

function arc_judge(arc::Arc,sign::Int64,bound::Bool,index::Int64)

    """
    arc_judge(arc::Arc,sign::Int64,bound::Bool,index::Int64)
    
    Check whether if the Line object is attached to the Arc vertex or not.

    # Arguments
    - `arc::Arc`: is the Arc object to be checked.
    - `sign::Int64`: is the parameter for interaction vertex (-1 is the left vertex & +1 is the right one).
    - `bound::Bool`: is the parameter determining the position of the Line relative to the vertex.(bound true is right, false is left)
    - `index::Int64`: is the Line index to be checked.
    """

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

    """
    extend!(diagram::Diagram)

    This function is not in use.
    Extend the last Line object in the diagram according to exponential distribution.

    # Arguments
    - `diagram::Diagram`: is the Diagram object to be modified.
    """

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

    """
    scale!(diagram::Diagram, order::Int64,m::Int64,μ::Float64,ω::Int64, samples::Int64)

    This function is not in use.
    Sample the diagram length according to Gamma distribution.

    # Arguments
    - `diagram::Diagram`: is the Diagram object to be modified.
    - `order::Int64`: is the number of Arc objects presented in the Diagram.
    - `m::Int64`: is the mass of the electron(assumed to be 1).
    - `μ::Float64`: is the energy offset for the Green's function(not in use).
    - `ω::Int64`: is the phonon frequency(assumed to be 1).
    - `samples::Int64`: is the number of times to sample the new diagram length.
    """

    line_box=diagram.line_box
    arc_box=diagram.arc_box
    end_arc_box=diagram.end_arc_box
    τ=diagram.τ
    record_τ=diagram.record_τ
    total_dis=diagram.dispersion
    
    coef=-τ/total_dis
    τ_new=rand(Gamma(2*order+1,coef),samples)

    if τ_new[1]>diagram.max_τ
        return false
    end

    diagram.record_τ=τ_new
    
    diagram.dispersion*=τ_new[1]/τ
    diagram.τ=τ_new[1]

    return true

end

function set_μ!(diagram::Diagram,new_μ::Float64)

    """
    set_μ!(diagram::Diagram,new_μ::Float64)

    This function is not in use.
    Set the new energy offset for the Diagram object.

    # Arguments
    - `diagram::Diagram`: is the Diagram object to be modified.
    - `new_μ::Float64`: is the new energy offset to be applied.
    """
    old_μ=copy(diagram.μ)
    total_dis=diagram.dispersion
    τ=diagram.τ

    diagram.dispersion+=old_μ*τ
    diagram.dispersion-=new_μ*τ
    diagram.μ=new_μ

    return diagram
end

function set_τ!(diagram::Diagram,new_τ::Float64)

    """
    set_τ!(diagram::Diagram,new_τ::Float64)

    Set the new diagram length for the Diagram object.(scaling with vertex proportion the same)

    # Arguments
    - `diagram::Diagram`: is the Diagram object to be modified.
    - `new_τ::Float64`: is the new length to be applied.
    """
    τ=diagram.τ
    diagram.dispersion/=τ
    diagram.dispersion*=new_τ
    diagram.τ=new_τ

    return diagram
end

function energy(diagram::Diagram)

    """
    energy(diagram::Diagram)

    Evaluate the ground energy for the Diagram object.

    # Arguments
    - `diagram::Diagram`: is the Diagram object to be evaluated.
    """
    total_dis=-diagram.dispersion
    τ=diagram.τ
    order=diagram.order
    return (total_dis-2*order)/τ
end

function mass_estimator(diagram::Diagram)

    """
    mass_estimator(diagram::Diagram)

    Evaluate the effective mass for the Diagram object.

    # Arguments
    - `diagram::Diagram`: is the Diagram object to be evaluated.
    """
    line_box=diagram.line_box
    p_squared = [0.0, 0.0, 0.0]
    τ=diagram.τ
    m=diagram.mass
    μ=diagram.μ
    for line in line_box
        p_squared .+= line.k*(line.period[2]-period[1])
    end
    return (norm(p_squared)^2)
end

function total_dis_check(diagram::Diagram)

    """
    total_dis_check(diagram::Diagram)

    Check the recorded total dispersion in update matches with the one calculated directly. (debug)

    # Arguments
    - `diagram::Diagram`: is the Diagram object to be evaluated.
    """

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
    
    """
    resample_arc!(diagram::Diagram,order::Int64,m::Int64,μ::Float64,ω::Int64,α_squared::Float64)

    This function is not in use.
    Resample one inserted Arc parameters (except the left vertex) according to the same steps in insert_arc! update.

    # Arguments
    - `diagram::Diagram`: is the Diagram object to be modified.
    - `order::Int64`: is the number of Arc objects presented in the Diagram.
    - `m::Int64`: is the mass of the electron(assumed to be 1).
    - `μ::Float64`: is the energy offset for the Green's function(not in use).
    - `ω::Int64`: is the phonon frequency(assumed to be 1).
    - `α_squared::Float64`: is 2pi*sqrt(2) of couopling strength.
    """
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
    index_out_r=arc.index_out
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

            deleteat!(line_box, index_out_r)
            sign_to_add=[sign_box[index_out_r-1][1],sign_box[index_out_r][2]]
            deleteat!(sign_box, index_out_r-1:index_out_r)
            insert!(sign_box, index_out_r-1, sign_to_add)
            if closed_arc
                if index_out_in>index_out_r
                    index_out_in-=1
                    new_arc.index_out=index_out_in
                    insert!(line_box, index_out_in, line_tem_3)
                    sign_to_add=[[sign_box[index_out_in-1][1],1],[1,sign_box[index_out_in-1][2]]]
                    deleteat!(sign_box,index_out_in-1)

                    for i in 1:2
                        insert!(sign_box, index_out_in-1, sign_to_add[3-i])
                    end
                elseif index_out_in<index_out_r && !cross_over
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
                    insert!(line_box, index_out_in, line_tem_3)
                    sign_to_add=[[sign_box[index_out_in-1][1],1],[1,sign_box[index_out_in-1][2]]]
                    deleteat!(sign_box,index_out_in-1)

                    for i in 1:2
                        insert!(sign_box, index_out_in-1, sign_to_add[3-i])
                    end
                end
            end
        end

        index_in=arc.index_in
        index_out=arc.index_out

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
            diagram.component-=1
        end

        #insert
        index_in=new_arc.index_in
        index_out=new_arc.index_out
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
            diagram.component+=1
        end

        return true
    end
end

function swap_arc!(diagram::Diagram)

    """
    swap_arc!(diagram::Diagram)

    Attempting to swap two nearest interaction vertexes of two Arcs in the Diagram object.

    # Arguments
    - `diagram::Diagram`: is the Diagram object to be modified.
    """
    order=diagram.order
    m=diagram.mass
    μ=diagram.μ
    ω=diagram.ω
    τ=diagram.τ

    # order checking
    if order<2
        return false
    end

    # randomly to choose Line to perform swap
    line_box=diagram.line_box
    line_index=rand(2:2*order)
    chosen_line=line_box[line_index]

    # covered property check
    if chosen_line.covered
        return false
    end

    # line length check to avoid abias in insert & remove (0.95 can be replaced with probability high enough e.g. 0.9-1)
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

    # find the Arc vertexes attached to the selected Line
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

    # debug to ensure covered property is set correctly
    if right_open && left_open && right_index == left_index
        println("overlap")
        return false
    end

    # Arcs to be swapped
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

    # new Arcs to be replaced after SWAP
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

    dis=(norm(chosen_line.k-sign[1]*q1+sign[2]*q2)^2-norm(chosen_line.k)^2)/(2m)
    dis+=(sign[1]-sign[2])*ω
    dis*=-(chosen_line.period[2]-chosen_line.period[1])*τ
    r=exp(dis) # acceptance ratio

    if r<rand()
        # return back to original diagram
        return false
    else
        # Arc SWAP is accepted. Update information in diagram
        diagram.dispersion+=dis

        # update line_box and sign_box
        deleteat!(line_box, line_index)
        insert!(line_box, line_index, new_line)
        sign_box=diagram.sign_box
        deleteat!(sign_box, line_index)
        insert!(sign_box, line_index, [sign[2],sign[1]])
        sign_box[line_index-1]=[sign_box[line_index-1][1],sign[2]]
        sign_box[line_index+1]=[sign[1],sign_box[line_index+1][2]]

        # arc_box update
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

        # covered property update
        if new_arc_l.index_out-new_arc_l.index_in == 2
            line_box[new_arc_l.index_in+1].covered=true
        end

        if new_arc_r.index_out-new_arc_r.index_in == 2
            line_box[new_arc_r.index_in+1].covered=true
        end
        sign=diagram.sign_box[line_index]

        return true
    end
end