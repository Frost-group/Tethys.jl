module Tethys

using StaticArrays # Gotta go quick

println("Loading Tethys...")

greet() = print("Hello Polarons!")

include("BareExpansion.jl")

greet()

end # module
