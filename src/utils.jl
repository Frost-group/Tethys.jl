include("measure.jl")
using CSV
using DataFrames
using JLD2,FileIO
using Logging

function slicematrix(A::AbstractMatrix)
    return [A[i, :] for i in 1:size(A,1)]
end

function get_data(directory, α, k)

    address=joinpath(directory,"α="*string(round(α, digits=3)),
            "k="*string(round(k, digits=3)))


    address = readdir(address; join=true)[1]  
    address = readdir(address; join=true)[1]

    df_total = DataFrame(CSV.File(joinpath(address,"total_green.csv"), header=0))
    df_zero = DataFrame(CSV.File(joinpath(address,"zero_green.csv"), header=0))

    diagram=load(joinpath(address,"diagram.jld2"), "diagram_a")
    hist=load(joinpath(address,"hist.jld2"), "hist_a")

    time_points=hist.time_points
    bin_width=hist.bin_width

    green_record = cumsum(slicematrix(Matrix(df_total)))
    zero_record = cumsum(slicematrix(Matrix(df_zero)))
    n_loop = length(green_record)

    normalized_data=hist.normalized_data
    bin_variance=jackknife(green_record,zero_record,n_loop,diagram,bin_width,0.1)


    return diagram, hist, green_record, zero_record, normalized_data, bin_variance
end

begin
    diagram, hist, green_record, zero_record, normalized_data, bin_variance = get_data("D://data",1.0,0.0 )
end
