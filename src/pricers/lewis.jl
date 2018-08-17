# Lewis, A. L. (2001)
# A simple option formula for general jump-diffusion and other exponential Lévy processes
# Available at SSRN 282110

struct Lewis <: AbstractPricer
  x::Array{Float64,1}
  w::Array{Float64,1}
  function Lewis(n=512)
    x, w = laguerre(n, 0.0)
    idx = w .!= 0.0
    new(x[idx], w[idx])
  end
end

_lewis_integrand(m, x0, x, k, r, q, T, u) =
  real(exp(u + logmgf(m, x0, T, 0.5 + u * im) + im * u * (x - k + (r - q) * T))) / (u^2 + 0.25)

function integral(m, p::Lewis, x0, x, k, r, q, T)
  z = 0.0
  @simd for i in 1:length(p.x)
    @inbounds z += p.w[i] * _lewis_integrand(m, x0, x, k, r, q, T, p.x[i])
  end
  z
end

function price(model, pricer::Lewis, x0, strike, rate, time; yield=0.0)
  stock = x0[1]
  stock * exp(-yield * time) - 1 / pi * sqrt(stock * strike) * exp(-0.5 * (rate + yield) * time) * integral(model, pricer, x0, log(stock), log(strike), rate, yield, time)
end
