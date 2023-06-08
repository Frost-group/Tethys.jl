include("measure.jl")
using CSV
using DataFrames
using JLD2,FileIO
using LsqFit
using Logging

function slicematrix(A::AbstractMatrix)
    return [A[i, :] for i in 1:size(A,1)]
end

function get_data(directory, α, k)

    address=joinpath(directory,"α="*string(round(α, digits=3)),
            "k="*string(round(k, digits=3)))

    diagram = []
    hist = []
    green_record = []
    zero_record = []
    normalized_data =[]
    bin_variance = []

    try
        address = readdir(address; join=true)[1]  
        address = readdir(address; join=true)[1]
        df_total = DataFrame(CSV.File(joinpath(address,"total_green.csv"), header=0))
        df_zero = DataFrame(CSV.File(joinpath(address,"zero_green.csv"), header=0))
    
        diagram=load(joinpath(address,"diagram.jld2"), "diagram_a")
        hist=load(joinpath(address,"hist.jld2"), "hist_a")

        bin_width=hist.bin_width
    
        green_record = cumsum(slicematrix(Matrix(df_total)))
        zero_record = cumsum(slicematrix(Matrix(df_zero)))
        n_loop = length(green_record)
    
        normalized_data=hist.normalized_data
        bin_variance=jackknife(green_record,zero_record,n_loop,diagram,bin_width,0.1)

    catch
        @info "Directory does not exist for "*address
    end



    return diagram, hist, green_record, zero_record, normalized_data, bin_variance
end

begin
    diagram, hist, green_record, zero_record, normalized_data, bin_variance = get_data("F://data",1.2,0.0 )
end


begin
    directory = "F://data"
    directory_list = readdir(directory; join=true)
    scanned_α = []
    energy_record = []
    E_error_record = []
    z0_record=[]
    Z_error_record=[]
    linear(t, p) = p[1].-p[2].*t

    for dir in directory_list
        α = parse(Float64, split(dir, "=")[2])
        @info "Loading α = "*string(α)
        diagram, hist, green_record, zero_record, green_func, bin_variance = get_data("F://data",α,0.0 )

        μ = diagram.μ
        bin_width = hist.bin_width
        min_time=Int(div(5,bin_width,RoundUp))
        max_time=Int(div(12,bin_width,RoundUp))

        time_points=hist.time_points[min_time:max_time]

        statis=sum(green_func[i,:] for i in 1:max_order+1)
        y=log.(statis)[min_time:max_time]
        w=1 ./(bin_variance[min_time:max_time]./(statis[min_time:max_time]).^2)

        p0=[0,(-α-1.26*(α/10)^2-μ)]
        fit = curve_fit(linear, time_points, y, w, p0)
        errors=standard_errors(fit)

        append!(scanned_α,α)
        append!(energy_record,fit.param[2]+μ)
        append!(E_error_record,errors[2])

        z0=exp(fit.param[1])
        append!(z0_record,z0)
        append!(Z_error_record,z0*errors[1])
    end

end
begin
    plot(scanned_α,energy_record,yerr=E_error_record,xlabel="α",ylabel="Energy",label="DiagMC")
    plot!(scanned_α,-scanned_α.-1.26*(scanned_α./10).^2,xlabel="α",ylabel="Energy",label="Pertub")
end

begin
    plot(scanned_α,z0_record,yerr=Z_error_record,xlabel="α",ylabel="Z_0",label="DiagMC")
end