import math

# --- Step 1: Define Input Parameters ---
# The problem provides the geometry and pavement type. We must make standard assumptions
# for parameters not explicitly given (Manning's n and rainfall intensity).

# Cross-slope (S) is given as 1.75%
S = 1.75 / 100

# The flow path length (L) is the width of the three lanes sloping in one direction.
# Lane width = 3.6 m, Number of lanes = 3.
L = 3.6 * 3

# Manning's roughness coefficient (n) for "rough-textured asphalt pavement".
# A standard value for this surface type is 0.016.
n = 0.016

# Rainfall intensity (i) in mm/hr. Since IDF curves are not provided, we must assume
# a reasonable value for a design storm that could lead to hydroplaning. A short-duration,
# high-intensity storm is appropriate. We assume an intensity of 150 mm/hr.
i = 150.0

# --- Step 2: Calculate Water Film Thickness ---
# We use the Gallaway formula, a standard empirical equation for pavement drainage design in SI units:
# Tw = 0.0038 * (n * L)^0.6 * i^0.4 / S^0.2

# Calculate each component of the formula
nl_val = n * L
nl_pow = math.pow(nl_val, 0.6)
i_pow = math.pow(i, 0.4)
s_pow = math.pow(S, 0.2)

# Calculate the final water film thickness (Tw)
Tw = 0.0038 * nl_pow * i_pow / s_pow

# --- Step 3: Display the Calculation and Result ---
# The final output shows the equation with the substituted values as requested.

print("Calculation of Design Water Film Thickness (Tw):")
print("Formula: Tw = 0.0038 * (n * L)^0.6 * i^0.4 / S^0.2")
print(f"Substituting the values:")
print(f"n (Manning's roughness) = {n}")
print(f"L (Flow path length) = {L:.1f} m")
print(f"i (Rainfall intensity) = {i:.1f} mm/hr")
print(f"S (Cross-slope) = {S}")
print("") # Blank line for spacing
print("Final Equation with numbers:")
print(f"Tw = 0.0038 * ({n} * {L:.1f})^0.6 * {i:.1f}^0.4 / {S}^0.2")
print(f"Tw = 0.0038 * ({nl_val:.4f})^0.6 * {i:.1f}^0.4 / {S:.4f}^0.2")
print(f"Tw = 0.0038 * {nl_pow:.4f} * {i_pow:.4f} / {s_pow:.4f}")
print(f"Tw = {(0.0038 * nl_pow * i_pow):.5f} / {s_pow:.5f}")
print(f"The design water film thickness is {Tw:.2f} mm.")