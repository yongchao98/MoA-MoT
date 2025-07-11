import numpy as np

# --- Step 1: Define constants and problem parameters ---
# Physical constants
c = 2.998e8  # Speed of light in m/s
a0 = 5.29e-11  # Bohr radius in m
hbar = 1.054e-34  # Reduced Planck constant in J*s
eps0 = 8.854e-12 # Permittivity of free space in F/m
e = 1.602e-19 # Elementary charge in C
alpha = e**2 / (4 * np.pi * eps0 * hbar * c) # Fine-structure constant

# Problem parameters
lambda_nm = 589  # Wavelength in nm
lambda_m = lambda_nm * 1e-9  # Wavelength in m
Z = 11  # Nuclear charge for Sodium
tau_exp_ns = 16.2  # Experimental lifetime in ns
tau_exp_s = tau_exp_ns * 1e-9  # Experimental lifetime in s

# From problem analysis
# Transition is 3p -> 3s, so l_max = max(1, 0) = 1
l_max = 1
# The hint g2/g1 = 2 points to the 3p_3/2 -> 3s_1/2 transition,
# for which the upper state degeneracy is g2 = 2*J+1 = 2*(3/2)+1 = 4.
g2 = 4

# --- Step 2: Calculate intermediate values ---
# Angular frequency
omega = 2 * np.pi * c / lambda_m

# Radial integral squared |I_r|^2
# From the provided wavefunctions, the integral evaluates to I_r = -27*sqrt(2)*a0/Z
# So, |I_r|^2 = (27**2 * 2) * (a0/Z)**2 = 1458 * (a0/Z)**2
Ir_sq = 1458 * (a0 / Z)**2

# --- Step 3: Calculate the theoretical transition rate A_total ---
# A_total = (4 * alpha * omega^3 * l_max * |I_r|^2) / (3 * c^2 * g2)
numerator = 4 * alpha * (omega**3) * l_max * Ir_sq
denominator = 3 * (c**2) * g2
A_total = numerator / denominator

# --- Step 4: Calculate the theoretical lifetime ---
tau_theo_s = 1 / A_total
tau_theo_ns = tau_theo_s * 1e9

# --- Step 5: Compare with experimental value ---
multiple = tau_theo_s / tau_exp_s

# --- Step 6: Print the results ---
print("--- Calculation of Theoretical Lifetime ---")
print(f"The formula for the theoretical lifetime is τ = 1 / A_total = (3 * c² * g₂) / (4 * α * ω³ * l_max * |I_r|²)")
print("\nPlugging in the values:")
print(f"c (speed of light) = {c:.3e} m/s")
print(f"g₂ (upper state degeneracy) = {g2}")
print(f"α (fine-structure constant) = {alpha:.4f}")
print(f"ω (angular frequency) = {omega:.3e} rad/s")
print(f"l_max (max orbital quantum number) = {l_max}")
print(f"|I_r|² (radial integral squared) = {Ir_sq:.3e} m²")

print(f"\nNumerator (3 * c² * g₂) = {3 * c**2 * g2:.3e}")
print(f"Denominator (4 * α * ω³ * l_max * |I_r|²) = {4 * alpha * omega**3 * l_max * Ir_sq:.3e}")

print(f"\nCalculated theoretical lifetime (τ_theo) = {tau_theo_s:.3e} s = {tau_theo_ns:.1f} ns")
print(f"Experimentally measured lifetime (τ_exp) = {tau_exp_s:.3e} s = {tau_exp_ns:.1f} ns")

print(f"\n--- Comparison ---")
print(f"The ratio of theoretical to experimental lifetime is {multiple:.2f}")
print("This means the calculated theoretical lifetime is approximately twice as long as the experimental value.")
<<<B>>>