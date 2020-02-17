# Starting with §3, a minimalist bare expansion for the Green function of the Frohlich polaron
# [1] Greitemann, J.; Pollet, L. Lecture Notes on Diagrammatic Monte Carlo for the Fröhlich Polaron. SciPost Phys. Lect. Notes 2018, 2. https://doi.org/10.21468/SciPostPhysLectNotes.2.

const MAX_ORDER=2
const α=1.5

struct Diag
    p # external momentum
    τ # external time

    O::Int # Order, integer

    phonon::SMatrix{MAX_ORDER, 3, Float64}
end


function BareExpansion(d::Diag)
# See (24) page 12 in [1], Trying to be as close to the maths as possible.
    
    # need to sort these in time in order to be able to calculate electron
    # momentum as phonons are emitted and absorbed
    # Generation of a phonon loses that phonon's momentum from the electron,
    # and vice versa

    # location of every phonon vertex in time, with change in momentum of
    # electron
    deltas=[ [p[1], -p[3], p[2], p[3]] for p in eachrow(d.phonon)]
    
    logG0tot=0
    p=d.p
    t=0
    for d in deltas
        Δτ=d[1]-t
        t=d[1]
        p=p+d[2]
        logG0tot+=logG0(p,Δτ)
    end

# God, what was I thinking of?
#for emit in d.phonon[]
#        deltas=emit[1],-emit[3]
#
#    emit=sortperm(d.phonon[1,:])
#    emitdeltas=[ [d.phonon[i,1],-d.phonon[i,3]] for i in emit]
#    
#    absorb=sortperm(d.phonon[2,:])
#    absorbdeltas=[ [d.phonons[i,2], +d.phonon[i,3]] for i in absorb]
#
#    deltas=sort(emitdeltas+absorbdeltas)


    GF=- α^d.O * 
    exp(logG0tot + sum(logD̃.(d.phonon[:,3], d.phonon[:,2] .- d.phonon[:,1])))
	# Is this all there is?
	
end

const μ=0 # chemical potential is the bane of my life
const m=1
function logG0(p,Δτ)
    -(p^2/2m-μ)*Δτ
end

const ωph=1
function logD̃(q,Δτ)
    -ωph * Δτ
end

### tests

#diag=Diag(0,1.0,MAX_ORDER,[0.1,0.2,1.0 , 0.3,0.4,1.0])

#println("Splish splash: ♓",diag,BareExpansion(diag))

