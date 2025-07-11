import math

# --- Given Parameters ---
# Band gap in eV
Eg = 3.0
# 1s exciton resonance peak energy in eV
E1s_peak = 1.0
# Principal quantum number for the target state
n = 3

# --- Step 1: Calculate the binding energy of the 1s exciton ---
# The exciton resonance peak is the band gap minus the binding energy.
# E_peak = Eg - E_binding  =>  E_binding = Eg - E_peak
Eb_1s = Eg - E1s_peak
print(f"Step 1: Calculate the 1s exciton binding energy.")
print(f"Binding Energy (1s) = Band Gap - Resonance Peak (1s)")
print(f"Eb_1s = {Eg} eV - {E1s_peak} eV = {Eb_1s} eV\n")


# --- Step 2: Calculate the effective 2D Rydberg energy (Ry_2D) ---
# For a 2D system, the binding energy is given by Eb_n = Ry_2D / (n - 0.5)^2
# For n=1, Eb_1s = Ry_2D / (1 - 0.5)^2 = Ry_2D / 0.25 = 4 * Ry_2D
# Therefore, Ry_2D = Eb_1s / 4
Ry_2D = Eb_1s / 4.0
print(f"Step 2: Calculate the effective 2D Rydberg energy.")
print(f"From Eb_1s = Ry_2D / (1 - 0.5)^2, we get Ry_2D = Eb_1s / 4.0")
print(f"Ry_2D = {Eb_1s} eV / 4.0 = {Ry_2D} eV\n")


# --- Step 3: Calculate the binding energy for the n=3 state ---
# "Rydberg energy for n=3" refers to the binding energy of the n=3 state.
# Eb_3s = Ry_2D / (3 - 0.5)^2
denominator = (n - 0.5)**2
Eb_3s = Ry_2D / denominator
print(f"Step 3: Calculate the Rydberg energy (binding energy) for n={n}.")
print(f"The final equation is: Eb_n = Ry_2D / (n - 0.5)^2")
print(f"Eb_{n} = {Ry_2D} / ({n} - 0.5)^2")
print(f"Eb_{n} = {Ry_2D} / {denominator}")
print(f"Result: Eb_{n} = {Eb_3s:.4f} eV")
