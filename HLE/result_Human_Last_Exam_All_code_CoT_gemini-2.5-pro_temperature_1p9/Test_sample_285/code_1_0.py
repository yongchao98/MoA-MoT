# The dimension of the space of coefficients a
n = 9

# The decay exponent of the integral I(a)
# We found this to be delta = 2/3
delta_numerator = 2
delta_denominator = 3
delta = delta_numerator / delta_denominator

print(f"The dimension of the space is n = {n}")
print(f"The decay exponent of the integral I(a) is delta = {delta_numerator}/{delta_denominator}")

# The integral of |I(a)|^p over R^n behaves like the integral of (r^(-delta))^p * r^(n-1) dr
# which is the integral of r^(n - 1 - p*delta) dr.
# This integral diverges if the exponent is >= -1.
# n - 1 - p*delta >= -1
# n >= p*delta
# p <= n / delta

p_numerator = n * delta_denominator
p_denominator = delta_numerator
p_critical = p_numerator / p_denominator

print(f"The condition for divergence is p <= n / delta.")
print(f"p <= {n} / ({delta_numerator}/{delta_denominator})")
print(f"p <= ({n} * {delta_denominator}) / {delta_numerator}")
print(f"p <= {p_numerator} / {p_denominator}")
print(f"p <= {p_critical}")

print(f"The largest value of p such that I is not in L^p(R^9) is {p_critical}.")