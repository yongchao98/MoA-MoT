import math

# The final formula for the ratio of differential cross-sections is R = 3 / (16 * cos^2(theta/2)).
# We are given the scattering angle theta.

# Given scattering angle
theta_rad = math.pi / 30

# The final formula involves three numerical parts: the numerator 3, the constant factor 16 in the denominator,
# and the term cos^2(theta/2). We will calculate and print each of these.

numerator = 3
denominator_const = 16
cos_squared_term = math.cos(theta_rad / 2)**2

# The total denominator is the product of the constant and the cosine term.
denominator_total = denominator_const * cos_squared_term

# The final ratio
ratio = numerator / denominator_total

print("Final equation: R = Numerator / (Denominator_Constant * cos^2(theta/2))")
print(f"Numerator: {numerator}")
print(f"Denominator_Constant: {denominator_const}")
print(f"theta (radians): {theta_rad}")
print(f"cos^2(theta/2): {cos_squared_term}")
print(f"Total Denominator: {denominator_total}")
print(f"The calculated ratio of differential cross-sections is: {ratio}")
