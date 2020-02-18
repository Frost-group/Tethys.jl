module Tethys

using StaticArrays # Gotta go quick
using Printf

println("Loading Tethys...")

greet() = println("Hello Polarons! ♓")

include("BareExpansion.jl")

greet()

end # module
