import math

# --- Step 1: Define constants and given values ---
q = 1.602e-19  # Electron charge in Coulombs
a0 = 5.292e-11  # Bohr radius in meters
c = 2.998e8   # Speed of light in m/s
epsilon0 = 8.854e-12 # Permittivity of free space in F/m
hbar = 1.0546e-34 # Reduced Planck constant in J.s
lambda_nm = 589  # Wavelength in nanometers
lambda_m = lambda_nm * 1e-9 # Wavelength in meters
tau_exp_ns = 16.2 # Experimental lifetime in nanoseconds
tau_exp_s = tau_exp_ns * 1e-9 # Experimental lifetime in seconds

# --- Transition details ---
# Transition is 3p -> 3s
n = 3
L_upper = 1  # 3p state
L_lower = 0  # 3s state
l_greater = max(L_upper, L_lower)
Z_eff = 1.0 # Effective nuclear charge for a valence electron

print("--- Calculating Theoretical Lifetime of Sodium 3p State ---")
print(f"Transition: {n}p -> {n}s")
print(f"Experimental lifetime (tau_exp) = {tau_exp_ns} ns")
print("\n--- Physical Constants and Parameters ---")
print(f"Wavelength (lambda) = {lambda_nm} nm")
print(f"Effective Nuclear Charge (Z_eff) = {Z_eff}")


# --- Step 2: Calculate intermediate quantities ---
# Angular frequency omega
omega = (2 * math.pi * c) / lambda_m

# Radial integral |<R_f|r|R_i>|
# Using the formula: (3*n / (2*Z)) * sqrt(n^2 - L^2) * a0
radial_integral = (3 * n / (2 * Z_eff)) * math.sqrt(n**2 - L_upper**2) * a0

# Line strength S_d = q^2 * l_greater * |radial_integral|^2
S_d = q**2 * l_greater * radial_integral**2

# Degeneracy of the upper orbital level (3p state)
g_L = 2 * L_upper + 1

print("\n--- Step-by-Step Calculation ---")
print(f"1. Angular Frequency (omega) = 2*pi*c / lambda = {omega:.4e} rad/s")
print(f"2. Radial Integral |<R_3s|r|R_3p>| = {radial_integral:.4e} m")
print(f"3. Line Strength (S_d) = q^2 * l_> * |radial_integral|^2 = {S_d:.4e} C^2 m^2")
print(f"4. Degeneracy of upper state (g_L = 2L+1) = {g_L}")


# --- Step 3: Calculate the Einstein A coefficient (spontaneous emission rate) ---
# A = (omega^3 * S_d) / (3 * pi * epsilon0 * hbar * c^3 * g_L)
A_numerator = omega**3 * S_d
A_denominator = 3 * math.pi * epsilon0 * hbar * c**3 * g_L
A_21 = A_numerator / A_denominator

print("\n--- Final Calculation of Lifetime ---")
print("Equation for emission rate A_21 = (omega^3 * S_d) / (3 * pi * epsilon0 * hbar * c^3 * g_L)")
print(f"Numerator = ({omega:.4e})^3 * {S_d:.4e} = {A_numerator:.4e}")
print(f"Denominator = 3 * pi * {epsilon0:.4e} * {hbar:.4e} * ({c:.4e})^3 * {g_L} = {A_denominator:.4e}")
print(f"Spontaneous Emission Rate (A_21) = {A_21:.4e} s^-1")


# --- Step 4: Calculate the theoretical lifetime ---
tau_th_s = 1 / A_21
tau_th_ns = tau_th_s * 1e9

print(f"Theoretical Lifetime (tau_th) = 1 / A_21 = {tau_th_s:.4e} s = {tau_th_ns:.2f} ns")


# --- Step 5: Compare theoretical and experimental lifetimes ---
ratio = tau_th_s / tau_exp_s
print(f"\n--- Comparison with Experiment ---")
print(f"The ratio of theoretical to experimental lifetime is: {tau_th_ns:.2f} ns / {tau_exp_ns} ns = {ratio:.2f}")

print("\nBased on this ratio, the calculated theoretical lifetime is approximately one-tenth (0.1) as long as the experimental measurement.")