logmgf(model, T, u) = ψ(model, u) * T - u * ψ(model, 1) * T

# BlackScholes
analyticinterval(::BlackScholes) = (-Inf, Inf)
ψ(m::BlackScholes, u) = 0.5 * m.σ^2 * u^2
