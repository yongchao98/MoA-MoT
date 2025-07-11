import math

# --- Given parameters ---
# Convert all units to SI (meters and kilograms)
m = 0.20  # kg, mass of the body
M = 0.80  # kg, mass of the guide
R = 0.20  # m (20 cm)
d = 0.50  # m (50 cm)
mu_D = 0.20  # friction coefficient

print("--- Step 1: Calculate the final height h using the work-energy theorem. ---")
# The formula is h = R - mu_D * d
h = R - mu_D * d
print(f"h = R - mu_D * d")
print(f"h = {R} m - {mu_D} * {d} m = {h:.2f} m")
print("-" * 50)

print("--- Step 2: Calculate the horizontal displacement of mass m relative to the guide. ---")
# The formula is Δx_m_rel = d + sqrt(2*h*R - h^2)
# The term under the square root is (2*h*R - h**2)
sqrt_term = math.sqrt(2 * h * R - h**2)
delta_x_m_rel = d + sqrt_term
print(f"Δx_m_rel = d + sqrt(2*h*R - h^2)")
print(f"Δx_m_rel = {d} + sqrt(2*{h:.2f}*{R} - {h:.2f}^2)")
print(f"Δx_m_rel = {d} + {sqrt_term:.4f} = {delta_x_m_rel:.4f} m")
print("-" * 50)

print("--- Step 3: Calculate the horizontal displacement of the guide. ---")
# The formula is ΔX_M = - (m / (m + M)) * Δx_m_rel
delta_X_M = - (m / (m + M)) * delta_x_m_rel

# Print the final calculation showing all the numbers
print("The displacement of the guide (ΔX_M) is calculated as:")
print(f"ΔX_M = - (m / (m + M)) * Δx_m_rel")
print(f"ΔX_M = - ({m} / ({m} + {M})) * {delta_x_m_rel:.4f}")
final_answer_m = delta_X_M
final_answer_cm = delta_X_M * 100

print("\n--- Final Result ---")
print(f"The horizontal displacement of the guide is: {final_answer_m:.4f} meters")
print(f"This is equivalent to {final_answer_cm:.2f} cm.")
print("The negative sign indicates that the guide moves to the left, opposite to the primary horizontal motion of the mass m.")

# Return the final answer in meters as requested by the format.
# Note: The prompt asks for the "displacement," which is a vector quantity, so the sign is important.
final_numeric_answer = final_answer_m
# The final answer in the format <<<...>>> should be a raw number
# Let's provide it with reasonable precision.
print(f"\n<<<d_final={final_numeric_answer}>>>")