# Starting with §3, a minimalist bare expansion for the Green function of the Frohlich polaron
# [1] Greitemann, J.; Pollet, L. Lecture Notes on Diagrammatic Monte Carlo for the Fröhlich Polaron. SciPost Phys. Lect. Notes 2018, 2. https://doi.org/10.21468/SciPostPhysLectNotes.2.
using Random

const MAX_ORDER=2
# See [1] 3.6 p.17 'Results'
const α=1.0
const μ=-1.2  # chemical potential is the bane of my life

struct Diag
    p # external momentum
    τ # external time

    O::Int # Order, integer

    phonon::MMatrix{MAX_ORDER, 3, Float64}
    move::MMatrix{MAX_ORDER, 3, Float64} # to hold move jumps
end


function BareExpansion(d::Diag; verbose=false)
# See (24) page 12 in [1], Trying to be as close to the maths as possible.
    
    # need to sort these in time in order to be able to calculate electron
    # momentum as phonons are emitted and absorbed
    # Generation of a phonon loses that phonon's momentum from the electron,
    # and vice versa

    # location of every phonon vertex in time, with change in momentum of
    # electron
    #deltas=[ (p[1], -p[3]),(p[2], p[3]) for p in eachrow(d.phonon)]
    
    # locate all el/ph interactions, and store the momentum exchange at this
    # point 
    deltas=vcat( [[p[1],p[3]] for p in eachrow(d.phonon)], [[p[2],-p[3]] for p in eachrow(d.phonon)])
    # time order the operators, so you can keep track of momentum
    timesorted=sortperm(deltas)

    logG0tot=0
    p=d.p
    t=0
    for i in timesorted
        ph=deltas[i]
        if verbose
            println("At time $(ph[1]), phonon momentum exchange $(ph[2])")
        end
        Δτ=ph[1]-t
        t=ph[1]
        p=p+ph[2]
        logG0tot+=logG0(p,Δτ)
    end
    # I don't know why, but this extra line cayses MC to diverge to -ve Inf
    # I'm certain there are off-by-one bugs in the above about line-up of sim
    # time with the phonon diagrams. But I am done for the day!
    #logG0tot+=logG0(p,d.τ-t) # last bare electron prop to end of sim time

    GF=- α^d.O * 
    exp(logG0tot + sum(logD̃.(d.phonon[:,3], d.phonon[:,2] .- d.phonon[:,1])))
	# Is this all there is?
end

const m=1
function logG0(p,Δτ)
#    println("Bare electron propagator, p= $(p) Δτ= $(Δτ)")
    -(p^2/2m-μ)*Δτ
end

const ωph=1
function logD̃(q,Δτ)
#    println("Bare Phonon propagator, q= $(q) Δτ= $(Δτ)")
    -ωph * Δτ
end

function Monte!(d::Diag;  verbose=false)
# Do you grow?
    GF0=BareExpansion(d)

    randn!(d.move)
    d.move ./= 100

    if verbose println("Moves: $(d.move)") end
    d.phonon .+= d.move

    # Brutal: checks to see whether *ANY* polaron value has gone negative
    # This is dumb as it also works on the polaron momentum
    if count(x->x<0, d.phonon)>0
        if verbose println(d.phonon, "Failed sanity test (time gone negative)") end
        d.phonon .-= d.move
        return GF0
    end

    GF1=BareExpansion(d)
    r=GF1/GF0

    if verbose println("GF0: $(GF0) GF1: $(GF1) r: $(r)") end

    if r > rand()
        if verbose println("Accept?!") end
        # log GW to histogram I guess...
        return GF1
    else
        if verbose println("Reject!") end
        d.phonon .-=  d.move  # delete our effects
        return GF0
    end

end

