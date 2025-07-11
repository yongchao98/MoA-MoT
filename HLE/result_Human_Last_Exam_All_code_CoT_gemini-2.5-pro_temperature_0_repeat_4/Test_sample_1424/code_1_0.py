import numpy as np

# Define the parameters
# omega: The frequency (energy) of the propagating particle (in eV)
omega = 2.5
# epsilon_k: The single-particle energy eigenvalue (in eV)
epsilon_k = 2.0
# eta: A small positive infinitesimal for causal boundary conditions
eta = 1e-9

# Calculate the denominator of the Green's function
denominator = (omega - epsilon_k) + 1j * eta

# Calculate the bare Green's function G_0
G0 = 1 / denominator

# Print the explanation and the final equation with all values
print("The bare Green's function G_0 is calculated using the formula:")
print("G_0(omega, epsilon_k) = 1 / (omega - epsilon_k + i*eta)\n")
print("For the given values:")
print(f"omega = {omega}")
print(f"epsilon_k = {epsilon_k}")
print(f"eta = {eta}\n")
print("The calculation is:")
# Using f-string formatting to display the full equation and result
# We format the complex number to show both real and imaginary parts clearly
print(f"G_0 = 1 / ({omega} - {epsilon_k} + i*{eta})")
print(f"G_0 = 1 / ({(omega - epsilon_k):.1f} + {eta}j)")
print(f"G_0 = {G0.real:.4f} + {G0.imag:.4f}j")