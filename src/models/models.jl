logmgf(model, x0, T, u) = ψ(model, u) * T - u * ψ(model, 1) * T

include("heston.jl")

# BlackScholes
analyticinterval(::BlackScholes) = (-Inf, Inf)
ψ(m::BlackScholes, u) = 0.5 * m.σ^2 * u^2

# Compound Poisson
struct CompoundPoisson
  λ::Float64
  μ::Float64
  δ::Float64
end
analyticinterval(m::CompoundPoisson) = (-3 / m.δ, 3 / m.δ)
ψ(m::CompoundPoisson, u) = m.λ * (exp(0.5 * m.δ^2 * u^2 + m.μ * u) - 1)

# Merton
ψ(m::Merton, u) =
  ψ(BlackScholes(m.α, m.σ), u) +
  ψ(CompoundPoisson(m.λ, m.μ, m.δ), u)
analyticinterval(m::Merton) = analyticinterval(CompoundPoisson(m.λ, m.μ, m.δ))

# Bates
function analyticinterval(m::Bates)
  a1, a2 = analyticinterval(Heston(m.α, m.κ, m.θ, m.σ, m.ρ))
  b1, b2 = analyticinterval(CompoundPoisson(m.λ, m.μ, m.δ))
  min(a1, b1), max(a2, b2)
end
logmgf(m::Bates, x0, T, u) =
  logmgf(Heston(m.α, m.κ, m.θ, m.σ, m.ρ), x0, T, u) +
  logmgf(CompoundPoisson(m.λ, m.μ, m.δ), x0, T, u)
