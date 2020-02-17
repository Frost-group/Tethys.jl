# Starting with §3, a minimalist bare expansion for the Green function of the Frohlich polaron
# [1] Greitemann, J.; Pollet, L. Lecture Notes on Diagrammatic Monte Carlo for the Fröhlich Polaron. SciPost Phys. Lect. Notes 2018, 2. https://doi.org/10.21468/SciPostPhysLectNotes.2.

const MAX_ORDER=2
const α=1.5

struct Diag
    p # external momentum
    τ # external time

    O::Int # Order, integer

    phonon::SMatrix{MAX_ORDER, 2, Float64}
end


function BareExpansion(d::Diag)
# See (24) page 12 in [1], Trying to be as close to the maths as possible.
    GF=- α^d.O
		
	
end


### tests

diag=Diag(0,1.0,5,[0.1,0.2 , 0.3,0.4])

println("Splish splash: ♓",diag,BareExpansion(diag))

