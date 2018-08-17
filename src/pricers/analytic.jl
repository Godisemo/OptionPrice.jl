struct Analytic <: AbstractPricer end

function price(model::BlackScholes, ::Analytic, x0, strike, rate, time; yield=0.0)
  stock = x0[1]
  d1 = (log(stock / strike) + (rate - yield + 0.5 * model.σ^2) * time) / (model.σ * sqrt(time))
  d2 = d1 - model.σ * sqrt(time)
  F = stock * exp((rate - yield) * time)
  exp(-rate * time) * (F * cdf(Normal(), d1) - strike * cdf(Normal(), d2))
end
