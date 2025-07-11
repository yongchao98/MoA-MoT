import math
from sympy import pi, Rational

# The sum is given by zeta(6) * zeta(8) / zeta(32).
# As zeta(32) is extremely close to 1, we can approximate the sum as zeta(6) * zeta(8).
# zeta(6) = pi^6 / 945
# zeta(8) = pi^8 / 9450

# The product is (pi^6 / 945) * (pi^8 / 9450)
numerator_zeta6 = 1
denominator_zeta6 = 945

numerator_zeta8 = 1
denominator_zeta8 = 9450

# Final expression is (numerator_zeta6 * numerator_zeta8) / (denominator_zeta6 * denominator_zeta8) * pi^(6+8)
final_numerator = numerator_zeta6 * numerator_zeta8
final_denominator = denominator_zeta6 * denominator_zeta8
final_pi_power = 14

# Simplify the fraction
common_divisor = math.gcd(final_numerator, final_denominator)
final_numerator //= common_divisor
final_denominator //= common_divisor

# The result is (1/8930250) * pi^14
print("The sum is given by the expression zeta(6) * zeta(8) / zeta(32).")
print("zeta(6) = pi^6 / 945")
print("zeta(8) = pi^8 / 9450")
print("zeta(32) is very close to 1, so we approximate the sum as zeta(6) * zeta(8).")
print(f"Sum approx zeta(6) * zeta(8) = (pi^6 / {denominator_zeta6}) * (pi^8 / {denominator_zeta8})")
print(f"= {final_numerator}/{final_denominator} * pi^{final_pi_power}")

# Final Answer
# To conform to the output format, print the full expression step-by-step
# Sum = ζ(6) * ζ(8)
# We can express these zeta functions in terms of pi.
# ζ(6) = (1/945) * π^6
# ζ(8) = (1/9450) * π^8
# So, the sum = (1/945) * π^6 * (1/9450) * π^8
# sum = (1 / (945 * 9450)) * π^14
# sum = (1 / 8930250) * π^14
print("Final Answer Calculation:")
print(f"{Rational(numerator_zeta6, denominator_zeta6)} * pi**6 * {Rational(numerator_zeta8, denominator_zeta8)} * pi**8")
print(f"= {Rational(final_numerator, final_denominator)} * pi**{final_pi_power}")
