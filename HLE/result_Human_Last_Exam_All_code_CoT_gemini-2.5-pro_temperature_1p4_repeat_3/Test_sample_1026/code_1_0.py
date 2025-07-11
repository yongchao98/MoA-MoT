import math

# Define the given variables
m = 0.20  # kg, mass of the body
M = 0.80  # kg, mass of the guide
R = 0.20  # m, radius of the circular arcs
d = 0.50  # m, length of the straight section
mu_D = 0.20 # coefficient of friction

# The core principle is the conservation of the horizontal position of the center of mass.
# The displacement of the guide (ΔX_M) can be related to the displacement of the
# mass m relative to the guide (Δx_rel) by the formula:
# ΔX_M = - (m / (m + M)) * Δx_rel
#
# The first half-oscillation takes the mass m from its horizontal starting position
# above one end of the straight section to its final horizontal position above the
# other end of the straight section.
# Therefore, the horizontal displacement of the mass relative to the guide is d.
delta_x_rel = d

# Calculate the horizontal displacement of the guide
# The negative sign indicates the guide moves in the opposite direction to the mass's relative motion.
guide_displacement = - (m * delta_x_rel) / (m + M)

# We can verify that the motion is possible by calculating the final height h.
# By conservation of energy: m*g*R = m*g*h + mu_D*m*g*d => h = R - mu_D*d
h = R - mu_D * d
if h > 0:
    print(f"The final height h is {h:.2f} m, which is positive, so the motion is possible.\n")
else:
    print("The mass does not reach the second arc, the problem premise is not met.\n")


# Print the final calculation step-by-step
print("The formula for the guide's displacement (ΔX_M) is:")
print("ΔX_M = - (m * d) / (m + M)\n")

print("Substituting the values:")
print(f"ΔX_M = - ({m} * {d}) / ({m} + {M})")

numerator = m * d
denominator = m + M
print(f"ΔX_M = - {numerator:.2f} / {denominator:.2f}\n")

print("Final Answer:")
print(f"The horizontal displacement of the guide is {guide_displacement:.2f} m.")

# Final answer in the required format
# <<<answer content>>>
# Note: The problem asks for the displacement. A displacement is a vector, so the sign is meaningful.
# It indicates the guide moves 0.10 m in the direction opposite to the mass's travel along the track.