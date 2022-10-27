using FFTW
using Logging
using Statistics

begin
    n_loop=1000
    α=1.5
    μ=-1.7
    num_mea=1; regime=Diff_more();
    p=0; max_τ=30; max_order=500; mass=1; ω=1;

    bin_width=max_τ/300
    min_time=Int(div(5,bin_width,RoundUp))
    max_time=Int(div(15,bin_width,RoundUp))

    hist=Hist_Record(300,max_τ,max_order)
    diagram=Diagram(p, max_τ, max_order, mass, μ, ω, α)
    green_record,green_func,variance=hist_measure_4!(diagram,hist,n_loop)
    @info "End of loop"

    time_points=hist.time_points[min_time:max_time]

    statis=sum(green_func[i,:] for i in 1:max_order+1)
    y=log.(statis)[min_time:max_time]

    F = fft(getindex.(green_record,1))
    acf = real.(ifft(F.*conj.(F)))
    acf = acf/(4*length(green_record))
    acf = acf/acf[1]
    #freqs = fftshift(fftfreq(length(t), fs))

    plot(collect(0: length(green_record)-1),acf)

    sum_threshold = Int(1000)
    ict = 1+2*sum(acf[1:sum_threshold])
    ref_ict = 5*ict

end