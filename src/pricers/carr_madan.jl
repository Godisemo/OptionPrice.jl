# Carr, P., & Madan, D. (1999)
# Option valuation using the fast Fourier transform
# Journal of computational finance, 2(4), 61-73

struct CarrMadan <: AbstractPricer
  x::Array{Float64,1}
  w::Array{Float64,1}
  function CarrMadan(n=512)
    x, w = laguerre(n, 0.0)
    idx = w .!= 0.0
    new(x[idx], w[idx])
  end
end

function _carr_madan_integrand(m, x, k, r, T, α, u)
  v = α + 1 + u * im
  real(exp(u + logmgf(m, T, v) + v * (x - k + r * T)) / (v^2 - v))
end

function integral(m, p::CarrMadan, x, k, r, T, α)
  z = 0.0
  @simd for i in 1:length(p.x)
    @inbounds z += p.w[i] * _carr_madan_integrand(m, x, k, r, T, α, p.x[i])
  end
  z
end

function price(model, pricer::CarrMadan, stock, strike, rate, time)
  logstock = log(stock)
  logstrike = log(strike)
  lb, ub = analyticinterval(model)
  f = α -> α * (rate * time + logstock - logstrike) + logmgf(model, time, α + 1) - log(α^2 + α)
  α = Optim.minimizer(optimize(f, max(0.0, lb - 1.0), min(1000.0, ub - 1.0)))
  exp(logstrike - rate * time) / pi * integral(model, pricer, logstock, logstrike, rate, time, α)
end
