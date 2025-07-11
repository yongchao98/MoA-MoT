import math

# This script presents the derived formula for the cutoff frequency at node a_0.

# The cutoff frequency (in radians per second) is given by omega_c = 1 / (R_th * C).
# Through analysis of the infinite ladder network, the Thevenin resistance R_th
# seen from node a_0 is found to be R_th = r * (1 + sqrt(3)).

# Substituting this into the formula for omega_c gives:
# omega_c = 1 / (r * C * (1 + sqrt(3)))

# We will now print the equation and the numerical constants within it,
# following the requested format.

numerator = 1
one_plus_sqrt_three_term_1 = 1
one_plus_sqrt_three_term_2 = 3

print("The formula for the angular cutoff frequency omega_c is of the form:")
print("omega_c = A / (r * C * (B + sqrt(D)))")
print("\nWhere the numerical constants are:")
print(f"A = {numerator}")
print(f"B = {one_plus_sqrt_three_term_1}")
print(f"D = {one_plus_sqrt_three_term_2}")

print("\nTherefore, the final equation is:")
# Using 'w_c' for omega_c
print(f"w_c = {numerator} / (r * C * ({one_plus_sqrt_three_term_1} + sqrt({one_plus_sqrt_three_term_2})))")
