import math

# 1. Define problem parameters.
# We can set the chord length c = 1.0, as the geometry is given in terms of c,
# and it will cancel out in the final ratio calculation.
c = 1.0
s = 0.5 * c
h = 0.5 * c

print(f"Given parameters (assuming chord c = {c}):")
print(f"Aerofoil separation, s = {s}")
print(f"Ride height, h = {h}")
print("-" * 30)

# 2. Define the coefficients K1 and K2 based on the derived aerodynamic interaction equations.
# K1 = 1 - c / (4 * h)
# K2 = c / (2 * s) - (c * s) / (2 * (s**2 + 4 * h**2))
# The lift ratio L1/L2 is given by the expression (K1 + K2) / (K1 - K2).

# 3. Calculate the value of K1.
k1_term_1 = 1.0
k1_term_2 = c / (4 * h)
K1 = k1_term_1 - k1_term_2

print("Calculating coefficient K1 = 1 - c / (4*h)")
print(f"K1 = {k1_term_1} - {c} / (4 * {h})")
print(f"K1 = {k1_term_1} - {k1_term_2}")
print(f"K1 = {K1}")
print("-" * 30)

# 4. Calculate the value of K2.
k2_term_1 = c / (2 * s)
k2_numerator = c * s
k2_denominator = 2 * (s**2 + 4 * h**2)
k2_term_2 = k2_numerator / k2_denominator
K2 = k2_term_1 - k2_term_2

print("Calculating coefficient K2 = c/(2*s) - (c*s) / (2*(s^2 + 4*h^2))")
print(f"First term of K2 = {c} / (2 * {s}) = {k2_term_1}")
print(f"Second term of K2 = ({c} * {s}) / (2 * ({s}**2 + 4 * {h}**2)) = {k2_numerator} / {k2_denominator:.4f} = {k2_term_2}")
print(f"K2 = {k2_term_1} - {k2_term_2}")
print(f"K2 = {K2}")
print("-" * 30)

# 5. Calculate the final lift ratio L1/L2.
ratio_numerator = K1 + K2
ratio_denominator = K1 - K2
lift_ratio = ratio_numerator / ratio_denominator

print("Calculating the final lift ratio L1/L2 = (K1 + K2) / (K1 - K2)")
print(f"L1/L2 = ({K1} + {K2}) / ({K1} - {K2})")
print(f"L1/L2 = {ratio_numerator} / {ratio_denominator}")
print(f"The final lift ratio L1/L2 is: {lift_ratio:.4f}")
