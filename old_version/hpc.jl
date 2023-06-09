include("Diagram.jl")
include("update.jl")
include("measure.jl")

using ArgParse

"""
    parse_commandline()

Parses input arguments from the command line to specify the simulation parameters.
Used for HPC and Archer2 job commands.

"""

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--alpha"
            help = "Coupling constant parameter"
            arg_type = Float64
            default = 1.0
        "--mu"
            help = "Chemical potential"
            arg_type = Float64
            default = -1.1
        "--k"
            help = "Electron momentum"
            arg_type = Float64
            default = 0.0
        "--sweeps"
            help = "Number of MC sweeps"
            arg_type = Int
            default = 10000
        "--updates"
            help = "Number of updates per loop (histograms sampled)"
            arg_type = Int
            default = 100000
        "--dir"
            help = "Directory for data storage"
            arg_type = String
            default = "D://data"
        "--save"
            help = "Flag to save data at specified directory"
            action = :store_true
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()

    α = parsed_args["alpha"]
    μ = parsed_args["mu"]
    p = parsed_args["k"]
    n_loop = parsed_args["sweeps"]
    directory = parsed_args["dir"]
    store_data = parsed_args["save"]

    max_τ=30; max_order=500; mass=1; ω=1;

    n_hist = 100000
    hist=Hist_Record(300,max_τ,max_order)
    diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)
    @info string(diagram.τ)
    diagram,hist,green_record,zero_record,green_func,variance=hist_measure!(diagram,hist,directory,n_loop,store_data,n_hist)

end

main()
