# Mass of the body in kg
m = 0.20
# Mass of the guide in kg
M = 0.80
# Radius of the circular arcs in cm
R_cm = 20
# Length of the straight section in cm
d_cm = 50

# Convert lengths from cm to meters for consistency in SI units
R_m = R_cm / 100
d_m = d_cm / 100

# The problem requires calculating the horizontal displacement of the guide (ΔX_M).
# The principle of conservation of the center of mass position is used.
# The formula derived is ΔX_M = - (m * (R + d)) / (m + M)
# The negative sign indicates that the guide moves in the opposite direction to the mass.

# Calculate the displacement of the guide
# Note: Δx_rel = R + d is the horizontal displacement of the mass relative to the guide
delta_x_rel = R_m + d_m
delta_X_M = - (m * delta_x_rel) / (m + M)

# Print the final equation with the numerical values
print("The formula for the guide's displacement (ΔX_M) is:")
print("ΔX_M = - (m * (R + d)) / (m + M)")
print("\nSubstituting the values (in SI units):")
print(f"ΔX_M = - ({m} * ({R_m} + {d_m})) / ({m} + {M})")
print(f"ΔX_M = - ({m} * {delta_x_rel}) / ({m + M})")
print(f"ΔX_M = - {m * delta_x_rel} / {m + M}")

# Print the final result
print(f"\nThe calculated horizontal displacement of the guide is {delta_X_M:.2f} meters.")
print("The negative sign indicates the displacement is to the left, opposite to the primary direction of the mass's motion.")

# The final answer in meters
final_answer = delta_X_M
# We can also output the absolute value as a magnitude, but displacement is a vector.
# Let's keep the signed value as it represents direction.
# The problem asks for the answer at the end.
print(f'<<<__{final_answer}__>>>')
