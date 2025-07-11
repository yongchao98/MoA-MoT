import numpy as np
import scipy.constants as const

# --- Step 1: Define constants and parameters (SI units) ---
# Provided in the problem
lambda_wl_nm = 589  # nm
tau_exp_ns = 16.2  # ns
Z = 11  # Nuclear charge of Sodium for the hydrogenic model

# Convert to SI units
lambda_wl = lambda_wl_nm * 1e-9  # m
tau_exp = tau_exp_ns * 1e-9  # s

# Physical constants from scipy
c = const.c
a0 = const.a0
alpha = const.alpha  # Fine-structure constant

print("--- Calculation of Theoretical Lifetime for Sodium 3p State ---")
print(f"Parameters: λ = {lambda_wl_nm} nm, Z = {Z}, Experimental τ = {tau_exp_ns} ns")
print("-" * 60)

# --- Step 2: Calculate the angular frequency of the transition ---
omega = 2 * np.pi * c / lambda_wl
print(f"Step 1: Angular Frequency (ω = 2πc/λ)")
print(f"ω = 2 * π * {c:.4e} / {lambda_wl:.4e} = {omega:.4e} rad/s\n")

# --- Step 3: Calculate the radial integral matrix element <3p|r|3s> ---
# The integral I = ∫ r³ * R_{3,1}(r) * R_{3,0}(r) dr can be solved analytically
# using the provided hydrogenic wavefunctions. The result is I = -3√2 * (a₀/Z).
radial_integral = -3 * np.sqrt(2) * (a0 / Z)
print(f"Step 2: Radial Integral (I = -3√2 * a₀/Z)")
print(f"I = -3 * √2 * {a0:.4e} / {Z} = {radial_integral:.4e} m\n")

# --- Step 4: Calculate the spontaneous emission rate (Einstein A coefficient) ---
# For 3p → 3s transition: l_upper = 1, l_lower = 0, l_max = 1.
l_upper = 1
l_max = 1
degeneracy_factor = l_max / (2 * l_upper + 1)

# Formula: A₂₁ = (4 * α * ω³) / (3 * c²) * (l_max / (2*l₂ + 1)) * I²
A21 = (4 * alpha * omega**3 / (3 * c**2)) * degeneracy_factor * radial_integral**2

print("Step 3: Spontaneous Emission Rate (A₂₁)")
print("A₂₁ = (4 * α * ω³ / (3 * c²)) * (l_max / (2*l₂+1)) * I²")
print(f"A₂₁ = (4 * {alpha:.4e} * ({omega:.4e})³ / (3 * ({c:.4e})²)) * ({l_max}/{2*l_upper+1}) * ({radial_integral:.4e})²")
print(f"A₂₁ = {A21:.4e} s⁻¹\n")

# --- Step 5: Calculate the theoretical lifetime ---
tau_th = 1 / A21
print("Step 4: Theoretical Lifetime (τ_th = 1/A₂₁)")
print(f"τ_th = 1 / {A21:.4e}")
print(f"τ_th = {tau_th:.4e} s  ({tau_th * 1e9:.1f} ns)\n")

# --- Step 6: Compare theoretical and experimental lifetimes ---
multiple = tau_th / tau_exp
print("Step 5: Comparison with Experimental Lifetime")
print(f"Ratio = τ_th / τ_exp = {tau_th * 1e9:.1f} ns / {tau_exp_ns} ns")
print(f"Ratio ≈ {multiple:.1f}\n")

print(f"The calculated theoretical lifetime is {multiple:.1f} times as long as the experimental one.")
print("The closest answer choice is 100.")
<<<E>>>