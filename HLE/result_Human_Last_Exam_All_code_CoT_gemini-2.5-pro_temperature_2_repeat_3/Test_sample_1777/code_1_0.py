import math

# Define the given physical constants
E_g = 3.0  # Band gap in eV
E_peak_1s = 1.0  # Resonance peak energy for the 1s exciton in eV
n = 3      # The principal quantum number of the target state

print("This script calculates the binding energy for an exciton in a specific state for a 2D material.")
print("-" * 30)

# Step 1: Calculate the binding energy of the 1s exciton (Eb_1s)
# The binding energy is the difference between the band gap and the exciton resonance peak energy.
# Formula: Eb = Eg - E_peak
Eb_1s = E_g - E_peak_1s
print("Step 1: Calculate the binding energy of the 1s (ground state) exciton (Eb_1s).")
print(f"Eb_1s = E_g - E_peak_1s = {E_g} eV - {E_peak_1s} eV = {Eb_1s} eV")
print("-" * 30)

# Step 2: Calculate the exciton Rydberg energy (Rx)
# For a 2D system, the binding energy formula is Eb_n = Rx / (n - 0.5)^2.
# For the ground state (n=1), this gives Eb_1s = Rx / (1 - 0.5)^2 = 4 * Rx.
# Therefore, Rx = Eb_1s / 4.
Rx = Eb_1s / 4.0
print("Step 2: Calculate the exciton Rydberg energy (Rx) using the 2D model.")
print(f"For the n=1 state, the formula Eb_1s = Rx / (1-0.5)^2 becomes Eb_1s = 4 * Rx.")
print(f"Therefore, Rx = Eb_1s / 4 = {Eb_1s} eV / 4 = {Rx} eV")
print("-" * 30)

# Step 3: Calculate the binding energy for the n=3 state (Eb_3)
# This is the "Rydberg energy for n=3" requested by the user.
# We use the 2D binding energy formula with n=3 and the calculated Rx.
Eb_3 = Rx / (n - 0.5)**2
print(f"Step 3: Calculate the binding energy for the n={n} state (Eb_{n}).")
print(f"Using the formula: Eb_{n} = Rx / (n - 0.5)^2")
print(f"Eb_{n} = {Rx} eV / ({n} - 0.5)^2 = {Rx} eV / {(n-0.5)**2} = {Eb_3:.4f} eV")

final_answer = Eb_3
print("\nThe Rydberg energy for n=3 is the binding energy of that state.")
print(f"Final Answer: {final_answer:.4f} eV")

print(f"\n<<<{final_answer:.4f}>>>")