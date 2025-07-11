import numpy as np

# This script calculates the force per unit area based on the formula from answer choice E,
# which is the most plausible answer despite discrepancies originating from apparent typos in the problem description.

# Define symbolic variables for the formula explanation
f_str = "f = (1/2) * (mu_0 * K_0^2 * cos^2(omega*t)) / (cosh^2(omega_p*d/c)) * exp(-omega*d/c)"

# Define numerical values for a sample calculation
mu_0 = 4 * np.pi * 1e-7  # Permeability of free space (H/m)
K_0 = 1000.0             # Surface current density amplitude (A/m)
omega = 1e9              # Angular frequency (rad/s), e.g., 1 GHz
t = 0.0                  # Time (s), choose t=0 for peak force
omega_p = 1e16           # Plasma frequency of the superconductor (rad/s)
d = 100e-9               # Thickness of the superconductor (m), e.g., 100 nm
c = 3e8                  # Speed of light in vacuum (m/s)

# --- Calculation step-by-step ---
val_one_half = 0.5
val_mu_0 = mu_0
val_K_0_sq = K_0**2
val_cos_sq = np.cos(omega * t)**2

# Argument of cosh
cosh_arg = omega_p * d / c
val_cosh_sq = np.cosh(cosh_arg)**2

# Argument of exp
exp_arg = -omega * d / c
val_exp = np.exp(exp_arg)

# Calculate final force magnitude
force_magnitude = val_one_half * (val_mu_0 * val_K_0_sq * val_cos_sq) / val_cosh_sq * val_exp

# --- Output the results ---
print("The formula for the force per unit area from choice E is:")
print(f"  {f_str}")
print("\nBelow are the values of each term in the equation for a sample scenario:")
print(f"  1/2 = {val_one_half}")
print(f"  mu_0 = {val_mu_0:.4e} H/m")
print(f"  K_0^2 = {K_0}^2 = {val_K_0_sq:.4e} A^2/m^2")
print(f"  cos^2(omega*t) = cos^2({omega:.1e}*{t}) = {val_cos_sq:.4f}")
print(f"  cosh^2(omega_p*d/c) = cosh^2({omega_p:.1e}*{d:.1e}/{c:.1e}) = {val_cosh_sq:.4e}")
print(f"  exp(-omega*d/c) = exp({exp_arg:.4f}) = {val_exp:.4f}")
print(f"\nFinal calculated force magnitude (at t={t} s): f_x = {force_magnitude:.4e} N/m^2")
