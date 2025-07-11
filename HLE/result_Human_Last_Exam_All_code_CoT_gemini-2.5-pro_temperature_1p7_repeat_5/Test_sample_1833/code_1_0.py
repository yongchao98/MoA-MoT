import math

# --- Given values ---
m = 1.0  # kg, mass of the ring
M = 1.0  # kg, mass of the object
g = 9.8  # m/s^2, acceleration due to gravity
theta_deg = 60.0  # degrees

# --- Calculations ---
# Convert angle to radians for use in trigonometric functions
theta_rad = math.radians(theta_deg)
sin_t = math.sin(theta_rad)
cos_t = math.cos(theta_rad)
cos_t_sq = cos_t**2

# Step 1: Calculate the L*omega^2 term from energy conservation
# L*omega^2 = 2*g*sin(theta)*(m+M) / (m + M*cos^2(theta))
l_omega_sq_numerator = 2 * g * sin_t * (m + M)
l_omega_sq_denominator = m + M * cos_t_sq
l_omega_sq = l_omega_sq_numerator / l_omega_sq_denominator

# Step 2: Calculate the tension T using Newton's second law
# T = (M*g*sin(theta) + M*L*omega^2) / (1 + (M/m)*cos^2(theta))
T_numerator_term1 = M * g * sin_t
T_numerator_term2 = M * l_omega_sq
T_numerator = T_numerator_term1 + T_numerator_term2
T_denominator = 1 + (M / m) * cos_t_sq
T = T_numerator / T_denominator

# --- Output the result step-by-step ---
print("Calculating the tension T in the string.")
print("The final formula for tension is: T = (M*g*sin(theta) + M*L*omega^2) / (1 + (M/m)*cos^2(theta))")
print("\nFirst, we calculate the term L*omega^2:")
print(f"L*omega^2 = (2 * {g} * sin({theta_deg}) * ({m} + {M})) / ({m} + {M} * cos^2({theta_deg}))")
print(f"L*omega^2 = ({l_omega_sq_numerator:.2f}) / ({l_omega_sq_denominator:.2f}) = {l_omega_sq:.2f}")

print("\nNext, we plug this into the tension formula:")
print(f"T = ({M}*{g}*sin({theta_deg}) + {M}*{l_omega_sq:.2f}) / (1 + ({M}/{m})*cos^2({theta_deg}))")
print(f"T = ({T_numerator_term1:.2f} + {T_numerator_term2:.2f}) / ({T_denominator:.2f})")
print(f"T = {T_numerator:.2f} / {T_denominator:.2f}")
print(f"\nThe calculated tension is: {T:.2f} N")

# Final answer in the required format
# print(f'<<<{T:.2f}>>>')