import math

# Define physical constants and parameters with example values
mu_0 = 4 * math.pi * 1e-7  # Permeability of free space (H/m)
R = 0.05                    # Radius of the outer coil (m)
N = 150                     # Turns in the second coil
N_0 = 600                   # Turns in the first coil
I_0 = 5.0                   # DC current (A)
i_0 = 0.5                   # Amplitude of AC current (A)
omega = 100 * math.pi       # Angular frequency (rad/s), for 50 Hz
g = 0.0015                  # Radial gap (m)
B_s = 1.8                   # Saturation flux density (T)
alpha_T = 2e-4              # Temperature coefficient of permeability (1/K)
T = 363.15                  # Operating temperature (K, 90 C)
T_0 = 293.15                # Reference temperature (K, 20 C)

# --- Calculation based on Answer B ---

# 1. Calculate the temperature-corrected permeability factor
temp_factor = 1 - alpha_T * (T - T_0)
mu_prime_0 = mu_0 * temp_factor

# 2. Calculate the saturation term in the denominator
saturation_term = 1 + (mu_0 * N_0 * I_0) / (g * B_s)

# 3. Calculate the numerator of the main fraction
# Note: In option B, the numerator uses mu_prime_0 for the overall permeability
# and the denominator uses mu_0 for the saturation part.
numerator_val = mu_prime_0 * N_0 * I_0 * i_0

# 4. Calculate the denominator of the main fraction
denominator_val = g**2 * saturation_term

# 5. Calculate the overall amplitude of the force
# The formula is f_x(t) = -2*pi*R*N * (fraction) * sin(omega*t)
amplitude = -2 * math.pi * R * N * (numerator_val / denominator_val)

# Output the components and the final equation
print("--- Equation Parameters ---")
print(f"Force Amplitude (coefficient of sin): {amplitude:.4f} N/m (assuming force per unit length)")
print(f"Angular Frequency (coefficient of t): {omega:.4f} rad/s")
print("\n--- Final Instantaneous Force Equation ---")
print(f"Based on choice B, the instantaneous force f_x(t) is:")
# Output the final equation with each number explicitly shown
print(f"f_x(t) = {amplitude:.4f} * sin({omega:.4f} * t)")