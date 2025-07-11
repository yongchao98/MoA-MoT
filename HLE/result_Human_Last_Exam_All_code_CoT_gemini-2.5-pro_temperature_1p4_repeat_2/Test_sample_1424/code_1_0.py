import numpy as np

# --- Parameters ---
# omega: Frequency (energy of the probe particle) in arbitrary energy units
omega = 2.0
# epsilon_k: Single-particle energy eigenvalue in the same energy units
epsilon_k = 1.5
# mu: Chemical potential in the same energy units
mu = 0.5
# delta: Small positive infinitesimal for causality. We use a small number for the calculation.
delta = 0.1

# --- Calculation ---
# The energy of the state relative to the chemical potential
xi_k = epsilon_k - mu

# Determine the sign of the imaginary term.
# The term is +1 for states above the Fermi level (unoccupied) and -1 for states below (occupied).
sign_xi_k = np.sign(xi_k)
# Handle the case where the state is exactly at the Fermi level
if sign_xi_k == 0:
    sign_xi_k = 1 # By convention, treat as unoccupied

# The complex denominator of the Green's function
denominator = (omega - xi_k) + 1j * delta * sign_xi_k

# The bare Green's function G_0
G0 = 1 / denominator

# --- Output the final equation with numbers ---
print("The functional dependence of the bare Green's function G_0 on the single-particle energy epsilon_k is given by:")
print("G_0(k, omega) = 1 / (omega - (epsilon_k - mu) + i*delta*sgn(epsilon_k - mu))")
print("\nFor the given values:")
print(f"omega = {omega}")
print(f"epsilon_k = {epsilon_k}")
print(f"mu = {mu}")
print(f"delta = {delta}\n")

# Reconstruct the equation string with the numerical values
imag_part_sign = "+" if sign_xi_k > 0 else "-"
print("The final equation with these values is:")
print(f"G_0 = 1 / ( {omega} - ({epsilon_k} - {mu}) {imag_part_sign} i * {delta} )")

# Print the result of the calculation
print(f"\nThe calculated bare Green's function is G_0 = {G0.real:.4f} {G0.imag:+.4f}j")