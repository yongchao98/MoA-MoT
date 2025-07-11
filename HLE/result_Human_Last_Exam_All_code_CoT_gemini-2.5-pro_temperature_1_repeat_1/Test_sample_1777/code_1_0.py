import math

# --- Given Parameters ---
# Band gap in eV
E_g = 3.0
# Resonance peak for the 1s exciton in eV
E_resonance_1s = 1.0
# Target principal quantum number
n = 3

# --- Calculation ---

# 1. Find the binding energy of the 1s state (n=1).
# Eb_1s = E_g - E_resonance_1s
Eb_1s = E_g - E_resonance_1s

# 2. Find the exciton Rydberg energy constant, Ry, from the 1s state.
# For a 2D exciton, Eb_1s = Ry / (1 - 0.5)^2 = Ry / 0.25 = 4 * Ry.
Ry = Eb_1s / 4.0

# 3. Calculate the binding energy for the n=3 state.
# This is the "Rydberg energy for n=3".
n_minus_half = n - 0.5
denominator = n_minus_half**2
Eb_n = Ry / denominator

# --- Output the Result ---
# The final output displays the full equation with the calculated numbers.

print(f"The Rydberg energy for n={n} is the binding energy of that state, Eb({n}).")
print(f"First, the exciton Rydberg constant (Ry) is found using the 1s state information:")
print(f"Ry = (E_gap - E_resonance_1s) / (1 - 0.5)^2 = ({E_g} - {E_resonance_1s}) / 0.25 = {Ry} eV")
print("\nThen, the binding energy for n=3 is calculated:")
print(f"Eb({n}) = Ry / (n - 0.5)^2")
print(f"Eb({n}) = {Ry} / ({n} - 0.5)^2 = {Ry} / {denominator} = {Eb_n:.3f} eV")
