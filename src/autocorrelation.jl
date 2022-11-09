include("Diagram.jl")
include("update.jl")
include("measure.jl")
using FFTW
using Logging
using Statistics

begin
    n_loop=1000
    n_hist = 10000
    cutoff = floor(Int, round(n_loop/10))
    α=1.5
    μ=-1.7
    num_mea=1; regime=Diff_more();
    p=0; max_τ=30; max_order=500; mass=1; ω=1;

    store_data = false

    bin_width=max_τ/300
    min_time=Int(div(5,bin_width,RoundUp))
    max_time=Int(div(15,bin_width,RoundUp))

    hist=Hist_Record(300,max_τ,max_order)
    diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)
    diagram,hist,green_record,zero_record,green_func,variance=hist_measure!(diagram,hist,"D://data",n_loop,store_data,n_hist)
    @info "End of loop"

    time_points=hist.time_points[min_time:max_time]

    statis=sum(green_func[i,:] for i in 1:max_order+1)
    y=log.(statis)[min_time:max_time]
    
    green_record_renorm = [green_record[i]./sum(green_record[i]) for i in 1:n_loop]

    F = fft(getindex.(green_record_renorm[cutoff:end],1))
    acf = real.(ifft(F.*conj.(F)))
    acf = acf/(4*length(green_record_renorm[cutoff:end]))
    acf = acf/acf[1]
    #freqs = fftshift(fftfreq(length(t), fs))

    display(plot(collect(0: length(green_record_renorm[cutoff:end])-1),acf))

    sum_threshold = Int(n_loop-cutoff)
    ict = 1+2*sum(acf[1:sum_threshold])
    ref_ict = 5*ict

end