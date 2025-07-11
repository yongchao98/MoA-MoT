import numpy as np

# This script calculates the bare Green's function G_0 for a given set of parameters.

# --- Theoretical Formula ---
print("In the Feynman path integral formalism, the bare Green's function G_0(k, ω)")
print("describes the propagation of a non-interacting particle.")
print("Its functional dependence on the single-particle energy eigenvalue ε_k is:\n")
print("  G_0(k, ω) = 1 / (ω - ε_k + iη * sgn(ε_k - μ))\n")
print("This shows that G_0 is inversely proportional to (ω - ε_k), with a small imaginary part for causality.\n")

# --- Example Calculation ---
print("Let's calculate G_0 for a specific case (an unoccupied 'particle' state):\n")

# Parameters
omega = 2.5      # Frequency (energy of the propagating particle) in eV
epsilon_k = 2.0  # Single-particle energy eigenvalue in eV
mu = 1.0         # Chemical potential (Fermi energy) in eV
eta = 1e-6       # A small positive infinitesimal

# --- Calculation Steps ---
# Determine the sign for the imaginary part
sign = np.sign(epsilon_k - mu)

# Calculate the denominator
denominator = (omega - epsilon_k) + 1j * eta * sign

# Calculate the Green's function
G0 = 1 / denominator

# --- Final Output ---
print("Using the following values:")
print(f"  ω = {omega}")
print(f"  ε_k = {epsilon_k}")
print(f"  μ = {mu}")
print(f"  η = {eta}\n")

print("The final equation with the numbers substituted is:")
print(f"G_0 = 1 / ({omega} - {epsilon_k} + 1j * {eta} * sgn({epsilon_k} - {mu}))")
print(f"G_0 = 1 / ({(omega - epsilon_k):.1f} + {1j * eta * sign})")
print(f"G_0 = {G0}")