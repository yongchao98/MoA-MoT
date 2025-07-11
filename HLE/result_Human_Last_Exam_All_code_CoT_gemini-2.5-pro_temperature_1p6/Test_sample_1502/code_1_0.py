import math

# This script demonstrates the conclusion of part (a).
# We showed that if p > 2*(1+3s)/(1+s), the functional J_t can become unbounded from below.
# Let's pick parameters that satisfy this condition and observe the behavior of J_t.

# --- Parameters ---
s = 2.0
# The critical value for p is 2*(1 + 3*s)/(1 + s)
p_critical = 2 * (1 + 3 * s) / (1 + s)
# We choose p > p_critical.
p = p_critical + 0.5  # So p = 14/3 + 0.5 ~ 5.17

# --- Functional Definition ---
# For demonstration, we simplify J_t by considering J(u_t, 0) and setting the
# values of the base integrals (for t=1) to 1.0.
# J_t(u_t, 0) = 0.5 * ||u_t||_{H^{1,s}}^2 - (mu_1/p) * ||u_t||_{Lp}^p
# We assume ||u||_{H^{1,s}}^2 = KineticPart + L2NormPart
# where KineticPart = integral(|d_x u|^2 + |d_y^s u|^2) and L2NormPart = integral(|u|^2)

kinetic_u = 1.0      # Placeholder value for the kinetic part of the H^{1,s} norm at t=1
l2_norm_sq_u = 1.0   # Placeholder for L2 norm squared (preserved under scaling)
lp_norm_p_u = 1.0    # Placeholder for Lp norm to the p at t=1
mu_1 = 2.0           # A positive coefficient

# --- Scaling Exponents ---
# Positive energy term (kinetic part) scales with t^(2s)
gamma_E = 2 * s
# Negative nonlinear term scales with t^gamma_p
gamma_p = (1 + s) * (p / 2 - 1)

# --- Coefficients for the J_t equation ---
# J(t) = C1 * t^gamma_E + C2 - C3 * t^gamma_p
term1_coeff = 0.5 * kinetic_u
term2_const = 0.5 * l2_norm_sq_u
term3_coeff = (mu_1 / p) * lp_norm_p_u

# --- Calculation and Output ---
print("Demonstrating unboundedness for J_t as t -> +inf.")
print(f"Chosen parameters: s = {s}, p = {p:.2f}")
print(f"Condition for unboundedness: p > {p_critical:.2f}, which is satisfied.")
print("\nFunctional form (derived from scaling laws):")
# Here we print the equation with the specific numbers calculated from our parameters
print(f"J(t) = {term1_coeff:.2f} * t**{gamma_E:.2f} + {term2_const:.2f} - {term3_coeff:.2f} * t**{gamma_p:.2f}")
print("-" * 40)

print("Value of J_t at different t:")
print("t \t\t J_t")
print("-" * 40)
for t_val in [1.0, 2.0, 4.0, 8.0, 16.0]:
    j_t = term1_coeff * (t_val**gamma_E) + term2_const - term3_coeff * (t_val**gamma_p)
    print(f"{t_val:<10.1f} \t {j_t:.3e}")

print("-" * 40)
print(f"As t increases, J_t becomes increasingly negative because the exponent of the negative term ({gamma_p:.2f})")
print(f"is larger than the exponent of the positive energy term ({gamma_E:.2f}).")
