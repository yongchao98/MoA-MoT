import math

# --- Given Parameters ---
m = 0.20      # Mass of the body in kg
M = 0.80      # Mass of the guide in kg
R = 0.20      # Radius of the circular arcs in meters (20 cm)
d = 0.50      # Length of the straight section in meters (50 cm)
mu_D = 0.20   # Coefficient of friction on the straight section

# --- Calculation ---
# The final formula for the guide's displacement is:
# ΔX_M = - (m / (m + M)) * Δx_rel
# where Δx_rel is the horizontal displacement of the mass relative to the guide.
# Δx_rel = d + sqrt(R^2 - (μ_D*d)^2)

# Step 1: Calculate the terms needed for the equation.
total_mass = m + M
mu_d_times_d = mu_D * d
sqrt_term = math.sqrt(R**2 - mu_d_times_d**2)
delta_x_rel = d + sqrt_term
displacement_guide = - (m / total_mass) * delta_x_rel

# --- Output ---
# As requested, printing the final equation with all numbers plugged in.
print("The final equation for the guide's displacement (ΔX_M) is:")
print(f"ΔX_M = - (m / (m + M)) * (d + sqrt(R² - (μ_D * d)²))")

print("\nSubstituting the given values:")
# Showing the equation with all the numbers
final_equation_str = f"ΔX_M = - ({m} / ({m} + {M})) * ({d} + sqrt({R}**2 - ({mu_D} * {d})**2))"
print(final_equation_str)

print("\nStep-by-step calculation:")
print(f"ΔX_M = - ({m} / {total_mass}) * ({d} + sqrt({R**2} - {mu_d_times_d**2}))")
print(f"ΔX_M = - {m / total_mass:.2f} * ({d} + {sqrt_term:.4f})")
print(f"ΔX_M = - {m / total_mass:.2f} * {delta_x_rel:.4f}")

# Final result
print("\nThe final calculated horizontal displacement of the guide is:")
print(f"{displacement_guide} meters")
