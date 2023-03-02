

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

    if right_open && left_open && right_index == left_index
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
    total_dis=dispersion(new_line, m, μ)-dispersion(chosen_line, m, μ)
    total_dis+=arc_dispersion(new_arc_r,ω)+arc_dispersion(new_arc_l,ω)-arc_dispersion(arc_r,ω)-arc_dispersion(arc_l,ω)
    r=exp(total_dis*diagram.τ)
    # println(w_y/w_x)
    # r=ratio*(1+log(ratio)/diagram.dispersion)^(2order+1)

    if r<rand()
        return false
    else
        diagram.dispersion+=total_dis*diagram.τ#log(r)
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
        update_arcp!(diagram,order,right_index,!right_open,m,μ)
        update_arcp!(diagram,order,left_index,!left_open,m,μ)
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
    diagram.dispersion-=dispersion(line, m, μ)*τ+dispersion(line_next, m, μ)*τ
    diagram.dispersion-=arc_dispersion(arc, ω)*τ
    τ_a=line.period[1]*τ
    τ_c=line_next.period[2]*τ
    e=(norm(line.k)^2-norm(line_next.k)^2)/2m+sign*ω
    τ_b=τ_a-log(1-rand()*(1-exp(-e*(τ_c-τ_a))))/e
    line.period[2]=τ_b/τ
    line_next.period[1]=τ_b/τ
    total_dis=0

    if sign==-1
        arc.period[1]=τ_b/τ
    else
        arc.period[2]=τ_b/τ
    end
    diagram.dispersion+=dispersion(line, m, μ)*τ+dispersion(line_next, m, μ)*τ
    diagram.dispersion+=arc_dispersion(arc, ω)*τ
end