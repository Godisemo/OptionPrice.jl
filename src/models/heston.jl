# Albrecher, H., Mayer, P., & Schoutens, W. T. J. (2007)
# The Little Heston Trap
# Wilmott Magazine, January issue, 83-92.

function analyticinterval(m::Heston)
  # an approximation of the analytic interval is the
  # area between the roots of p(x) = ax^2+bx+c
  a = (1.0 - m.ρ^2) * m.σ^2
  b = 2.0 * m.ρ * m.σ * m.κ - m.σ^2
  c = -m.κ^2
  A = -0.5 * b / a
  B = sqrt(0.25 * b^2 / a^2 - c / a)
  A - B, A + B
end

function logmgf(m::Heston, x0, T, u)
  v0 = x0[2]
  σ2 = m.σ * m.σ
  a0 = m.ρ * m.σ * u - m.κ
  a1 = sqrt(a0 * a0 + σ2 * u * (1 - u))
  a2 = m.ρ * m.σ * u - m.κ
  a3 = -(a2 + a1)
  a4 = a3 / (a1 - a2)
  a5 = exp(-a1 * T)
  a6 = 1 - a4 * a5
  ((a3 * T - 2 * log(a6 / (1 - a4))) * m.θ * m.κ + v0 * a3 * (1 - a5) / a6) / σ2
end
