import cmath
import numpy as np

# In the Feynman path integral formalism at finite temperature, physical quantities
# are often calculated in imaginary time and then transformed to Matsubara frequency space.

# --- Parameters for our example ---

# A single-particle energy eigenvalue (relative to the chemical potential).
# This represents the energy of a quantum state 'k'.
epsilon_k = 2.0

# In the imaginary time formalism, we use Matsubara frequencies.
# For fermions, these are given by omega_n = (2n+1) * pi / beta,
# where n is an integer and beta is the inverse temperature (1/k_B T).
# Let's choose the simplest case: n=0 and beta=1.0 for simplicity.
n = 0
beta = 1.0
omega_n = (2 * n + 1) * np.pi / beta

# --- The Functional Dependence ---

# The bare Green's function G_0(k, i*omega_n) in Matsubara frequency space
# has a fundamental dependence on the single-particle energy epsilon_k.
# The formula is: G_0 = 1 / (i*omega_n - epsilon_k)
# This shows that G_0 is inversely proportional to the difference between
# the complex frequency (i*omega_n) and the single-particle energy (epsilon_k).

# --- Calculation ---

# Calculate the complex frequency term
i_omega_n = 1j * omega_n

# Calculate the bare Green's function using the formula
G0 = 1 / (i_omega_n - epsilon_k)

# --- Output ---

print("The bare Green's function G_0 describes a non-interacting particle.")
print("In frequency space, its functional dependence on the single-particle energy eigenvalue epsilon_k is:\n")
print("G_0(k, i*omega_n) = 1 / (i*omega_n - epsilon_k)\n")
print("This shows G_0 is inversely proportional to the term (i*omega_n - epsilon_k).")

print("\n--- Example Calculation ---")
print(f"Let's use a single-particle energy epsilon_k = {epsilon_k}")
print(f"And the lowest fermionic Matsubara frequency i*omega_n (for n=0, beta=1) = {i_omega_n:.4f}")
print("\nThe final equation with these numbers is:")
# The following line explicitly prints the equation with the values, as requested.
print(f"G_0 = 1 / ({i_omega_n:.4f} - {epsilon_k:.4f})")
print(f"\nThe resulting bare Green's function is: G_0 = {G0:.4f}")
