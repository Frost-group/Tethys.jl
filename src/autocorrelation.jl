using FFTW
using Logging
using Statistics

"""
  auto_window(taus, c)

Automated window scheme proposed by Sokal, where c=5.
See A. Sokal, Monte Carlo Methods in Statistical Mechanics: Foundations and New Algorithms. 1997.

"""
function auto_window(taus, c::Float64)
    m = collect(1:length(taus)) .< c * taus
    if sum(m) > 0
        return argmin(m)
    end
    return length(taus)
end

"""
  autocorrelation(x, burnin)

Calculates the autocorrelation of the 1D estimator x, given the burnin samples.
Computed using Fast Fourier Transforms (FFTs).

"""
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
