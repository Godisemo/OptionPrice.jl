using OptionPrice
using SDEModels
@static if VERSION < v"0.7.0-DEV.2005"
  using Base.Test
else
  using Test
end

# write your own tests here

analytic = Analytic()
lewis = Lewis()
carr_madan = CarrMadan()

@test price(BlackScholes(0.01, 0.1), analytic, 100.0, 90.0, 0.01, 2.0) ≈ 13.147550110559813
@test price(BlackScholes(0.01, 0.1), carr_madan, 100.0, 90.0, 0.01, 2.0) ≈ 13.147550110559813
@test price(BlackScholes(0.01, 0.1), lewis, 100.0, 90.0, 0.01, 2.0) ≈ 13.147550110559813

@test price(Heston(0.01, 10, 0.16, 0.1, -0.8), carr_madan, [90.0, 0.16], 80.0, 0.01, 1.0) ≈ 19.489631391639065
@test price(Heston(0.01, 10, 0.16, 0.1, -0.8), lewis, [90.0, 0.16], 80.0, 0.01, 1.0) ≈ 19.489631391639065
