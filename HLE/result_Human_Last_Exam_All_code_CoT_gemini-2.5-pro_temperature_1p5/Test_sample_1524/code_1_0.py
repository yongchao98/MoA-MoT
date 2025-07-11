import numpy as np

def bose_einstein_occupation(energy, mu, kBT):
  """
  Calculates the average Bose-Einstein occupation number for a given state.
  The equation is: n(ε) = 1 / (exp((ε - μ) / k_B*T) - 1)
  """
  exponent = (energy - mu) / kBT
  # If the argument to exp() is zero or negative, the result is unphysical
  # or divergent. This corresponds to μ >= ε.
  if exponent <= 1e-9: # Use a small epsilon to avoid float precision issues
    return float('inf')
  return 1.0 / (np.exp(exponent) - 1.0)

# --- Parameters for Demonstration ---
# For simplicity, we can set the ground state energy to 0.
# The choice of units is arbitrary; we can set k_B*T = 1.
ground_state_energy = 0.0 # This is ε_0
thermal_energy = 1.0      # This is k_B*T

print("--- Demonstrating the Limit on Chemical Potential (μ) for BEC ---")
print(f"We analyze the Bose-Einstein occupation number n(ε) for the ground state (ε_0).")
print(f"We will use the following values in the equation n(ε_0) = 1 / (exp((ε_0 - μ) / k_B*T) - 1):")
print(f"Ground State Energy (ε_0): {ground_state_energy}")
print(f"Thermal Energy (k_B*T):    {thermal_energy}\n")

print("The fundamental condition for a physical occupation number is that μ must be less than ε.")
print("For the ground state, this means μ < ε_0.")
print("Bose-Einstein Condensation occurs when the ground state becomes macroscopically occupied.")
print("This happens as the chemical potential μ approaches the ground state energy ε_0 from below.")
print("\nLet's observe the ground state occupation n(ε_0) as μ gets closer to ε_0:\n")

# A list of chemical potential values approaching the ground state energy
mu_values = [
    ground_state_energy - 0.1,
    ground_state_energy - 0.01,
    ground_state_energy - 0.001,
    ground_state_energy - 0.0001,
]

print("--------------------------------------------------")
print("|   μ (Chemical Potential)   |  n(ε_0) (Occupation) |")
print("--------------------------------------------------")
for mu_val in mu_values:
  # Calculate each number in the final equation
  energy = ground_state_energy
  mu = mu_val
  kBT = thermal_energy
  occupation = bose_einstein_occupation(energy, mu, kBT)
  print(f"| {mu:^26.4f} | {occupation:^20.2f} |")
print("--------------------------------------------------\n")

print("Conclusion:")
print("As μ approaches ε_0, the number of particles in the ground state n(ε_0) grows without bound.")
print("This allows for a macroscopic fraction of particles to condense into the ground state.")
print(f"Thus, the fundamental limit on the chemical potential is μ = ε_0.")
print("This value is precisely the chemical potential of a non-interacting Bose gas at T=0, as stated in option C.")