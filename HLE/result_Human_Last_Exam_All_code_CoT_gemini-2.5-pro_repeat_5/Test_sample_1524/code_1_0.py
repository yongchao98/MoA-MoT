import numpy as np

def bose_einstein_occupation(epsilon, mu, k_B_T=1.0):
    """Calculates the occupation number for a given energy state."""
    denominator = np.exp((epsilon - mu) / k_B_T) - 1
    if denominator <= 0:
        return float('inf') # Represents unphysical or infinite occupation
    return 1 / denominator

# --- Setup for Demonstration ---
# Let's define a ground state energy. The value is arbitrary.
epsilon_0 = 0.5  # Ground state energy in some units
# We'll use thermal energy k_B_T = 1.0 for simplicity.

print(f"Illustrating the limit on chemical potential (μ) for a system with ground state energy ε₀ = {epsilon_0}")
print("The occupation number n₀ = 1 / (exp((ε₀ - μ) / k_B*T) - 1)")
print("-" * 60)

# Case 1: μ is significantly less than ε₀
mu_1 = 0.1
n_1 = bose_einstein_occupation(epsilon_0, mu_1)
print(f"Equation with μ = {mu_1}: n₀ = 1 / (exp(({epsilon_0} - {mu_1}) / 1.0) - 1) = {n_1:.4f}")
print("Result: For μ far below ε₀, the ground state occupation is small.")

# Case 2: μ gets closer to ε₀
mu_2 = 0.49
n_2 = bose_einstein_occupation(epsilon_0, mu_2)
print(f"\nEquation with μ = {mu_2}: n₀ = 1 / (exp(({epsilon_0} - {mu_2}) / 1.0) - 1) = {n_2:.4f}")
print("Result: As μ approaches ε₀, the occupation grows rapidly.")

# Case 3: μ is very close to ε₀, simulating condensation
mu_3 = 0.4999
n_3 = bose_einstein_occupation(epsilon_0, mu_3)
print(f"\nEquation with μ = {mu_3}: n₀ = 1 / (exp(({epsilon_0} - {mu_3}) / 1.0) - 1) = {n_3:.4f}")
print("Result: For μ infinitesimally close to ε₀, the occupation becomes macroscopic (BEC).")

# Case 4: μ ≥ ε₀, which is unphysical
mu_4 = 0.5
n_4 = bose_einstein_occupation(epsilon_0, mu_4)
print(f"\nEquation with μ = {mu_4}: n₀ = 1 / (exp(({epsilon_0} - {mu_4}) / 1.0) - 1) = {n_4}")
print("Result: If μ = ε₀, the occupation is infinite (unphysical denominator is zero).")

print("-" * 60)
print("Conclusion: The chemical potential μ must be less than the ground state energy ε₀.")
print("The chemical potential of a non-interacting Bose gas at T=0 is exactly ε₀, which serves as the fundamental limit.")
