using FFTW
using Logging
using Statistics


function auto_window(taus, c::Float64)
    m = collect(1:length(taus)) .< c * taus
    if sum(m) > 0
        return argmin(m)
    end
    return length(taus)
end

function autocorrelation(x, burnin::Int64)

    n = Int64(log2(nextpow(2,length(x[burnin:end]))))

    F = fft(x[burnin:end].-mean(x[burnin:end]))
    acf = real.(ifft(F.*conj.(F)))
    acf = acf/(4*n)
    acf = acf/acf[1]

    taus = 2.0 * cumsum(acf) .- 1.0

    window = auto_window(taus, 5.0)
    return taus[window]
end
