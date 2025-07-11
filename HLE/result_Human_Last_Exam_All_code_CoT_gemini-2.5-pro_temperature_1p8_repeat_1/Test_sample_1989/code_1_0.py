import numpy as np

# User-defined parameters for demonstration
A = 1.0
B = 2.0
r = 100.0
theta_deg = 45.0  # Angle in degrees

# Convert theta to radians for calculations
theta_rad = np.deg2rad(theta_deg)

# --- Calculation of the Asymptotic Formula Components ---
# The large-distance behavior is of the form:
# C * r^(p) * exp(q) * exp(-r * (1 - cos(theta)))
# where p and q depend on A, B, and theta.

# 1. Base term (from Oseen solution, A=0, B=0)
p_base = -0.5

# 2. Corrector terms due to non-zero A and B
# The power of r has a correction:
p_corr = A * (1 - np.cos(theta_rad)) + B * np.sin(theta_rad)
# There is a new exponential (phase) term:
q_corr = B * theta_rad

# 3. Full power and phase terms
p_total = p_base + p_corr
q_total = q_corr
exp_decay_term = -r * (1 - np.cos(theta_rad))

# --- Output the results ---
print(f"Analysis of the asymptotic behavior for omega(r, theta)")
print("-" * 60)
print(f"For parameters: A = {A}, B = {B}, r = {r}, theta = {theta_deg} degrees")
print("-" * 60)
print("The asymptotic formula is proportional to:")
print("r^(p) * exp(q) * exp(-r*(1-cos(theta)))\n")
print("Here are the computed values for the exponents 'p' and 'q':")

print("\n1. Final power of r, 'p':")
print(f"p = -0.5 + A*(1 - cos(theta)) + B*sin(theta)")
print(f"p = -0.5 + {A:.1f}*(1 - {np.cos(theta_rad):.4f}) + {B:.1f}*{np.sin(theta_rad):.4f}")
print(f"p = {p_total:.4f}")

print("\n2. Final phase term exponent, 'q':")
print(f"q = B*theta")
print(f"q = {B:.1f} * {theta_rad:.4f}")
print(f"q = {q_total:.4f}")

print("\n3. Dominant exponential decay term:")
print(f"Value of exponent = -r*(1-cos(theta))")
print(f"Value = -{r:.1f}*(1-{np.cos(theta_rad):.4f}) = {exp_decay_term:.4f}")

print("\n--- Summary ---")
print("The large-distance behavior is proportional to:")
print(f"r^({p_total:.4f}) * exp({q_total:.4f}) * exp({exp_decay_term:.4f})")
