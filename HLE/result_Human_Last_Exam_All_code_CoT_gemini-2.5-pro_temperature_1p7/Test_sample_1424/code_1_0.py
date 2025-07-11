import cmath
import numpy as np

def calculate_bare_greens_function(omega, epsilon_k, mu=0.0, delta=1e-6):
    """
    Calculates the bare Green's function G_0(k, omega).

    Args:
        omega (float): The frequency of the particle.
        epsilon_k (float): The single-particle energy eigenvalue.
        mu (float): The chemical potential.
        delta (float): A small positive infinitesimal.

    Returns:
        complex: The value of the bare Green's function.
    """
    # Determine the sign of the imaginary part based on the position relative to the chemical potential
    sign = np.sign(epsilon_k - mu)
    
    # Handle the case where epsilon_k is exactly at the chemical potential
    if sign == 0:
        sign = 1 # Conventionally, states at the chemical potential are treated as unoccupied

    # The denominator of the Green's function
    denominator = (omega - epsilon_k) + 1j * delta * sign

    # The Green's function is the inverse of the denominator
    G0 = 1 / denominator

    return G0, sign

# --- Example Calculation ---
# You can change these values to see how the result changes.
# This example considers a particle in an unoccupied state (epsilon_k > mu).
omega_val = 2.5  # Energy/frequency of the propagator
epsilon_k_val = 2.0  # Energy of the single-particle state
mu_val = 1.0     # Chemical potential (Fermi level)
delta_val = 1e-6 # Small infinitesimal for causality

# Perform the calculation
G0_result, sign_val = calculate_bare_greens_function(omega_val, epsilon_k_val, mu_val, delta_val)

# --- Print the results ---
print("--- Bare Green's Function Calculation ---")
print(f"Inputs:")
print(f"  Particle Frequency (ω): {omega_val}")
print(f"  Single-Particle Energy (ϵ_k): {epsilon_k_val}")
print(f"  Chemical Potential (μ): {mu_val}")
print("-" * 20)

print("Formula:")
if sign_val > 0:
    print("  G_0 = 1 / (ω - ϵ_k + iδ)")
    sign_str = '+'
else:
    print("  G_0 = 1 / (ω - ϵ_k - iδ)")
    sign_str = '-'

print("\nCalculation with substituted values:")
# Print the final equation with the numbers filled in
imaginary_part_str = f"i*{delta_val}"
if sign_val < 0:
    imaginary_part_str = f"-i*{delta_val}"
    
print(f"  G_0 = 1 / ({omega_val} - {epsilon_k_val} {sign_str} {imaginary_part_str})")
print(f"  G_0 = 1 / ({(omega_val - epsilon_k_val)} + {1j*delta_val*sign_val})")


print("\nFinal Result:")
print(f"  G_0(k, ω) = {G0_result}")
