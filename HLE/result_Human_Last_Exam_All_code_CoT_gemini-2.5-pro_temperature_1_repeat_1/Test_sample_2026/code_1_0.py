import math

# The problem reduces to finding the surface area of a sphere defined by the equation:
# x0^2 + y0^2 + z0^2 = R^2

# From the problem derivation, the constant part of the sum (alpha + beta + gamma) is:
sum_const = 10**25

# The radius squared (R^2) of the sphere is given by:
# R^2 = sum_const / 2
R_squared = sum_const / 2

# The surface area of a sphere is A = 4 * pi * R^2
# Substituting R^2, we get A = 4 * pi * (sum_const / 2) = 2 * pi * sum_const
two = 2
four = 4
pi_val = math.pi
area = two * pi_val * sum_const

print("The solvability condition for the initial values (x0, y0, z0) results in the equation of a sphere:")
print(f"x0^2 + y0^2 + z0^2 = R^2")
print(f"where R^2 = {sum_const} / 2 = {R_squared:.1e}\n")

print("The area of this sphere is calculated using the formula A = 4 * pi * R^2.")
print("The final equation for the area, after substitution, is A = 2 * pi * 10^25.\n")
print("Here are the numbers in the final equation:")
print(f"Constant factor: {two}")
print(f"Value of pi: {pi_val}")
print(f"Power of 10 factor: {sum_const:.0e}\n")

print(f"The calculated area is: {area:.6e}")