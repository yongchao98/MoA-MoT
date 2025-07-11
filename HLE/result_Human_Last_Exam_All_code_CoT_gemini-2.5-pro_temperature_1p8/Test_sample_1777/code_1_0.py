import math

# --- Given parameters ---
# Band gap in eV
Eg = 3.0
# 1s exciton resonance peak in eV
E_ex_1s = 1.0
# Principal quantum number for the target state
n = 3

# --- Step 1: Calculate the 1s exciton binding energy (the exciton Rydberg energy) ---
# The relationship is E_ex = Eg - Eb, so Eb = Eg - E_ex.
Eb_1s = Eg - E_ex_1s

print("Step 1: Calculate the 1s exciton binding energy, Eb(1s).")
print(f"This is found by subtracting the 1s exciton peak energy from the band gap.")
print(f"Eb(1s) = Eg - E_ex(1s) = {Eg} eV - {E_ex_1s} eV = {Eb_1s} eV")
print("-" * 30)

# --- Step 2: Use the 2D exciton formula to find the binding energy for n=3 ---
# For a 2D system, the binding energy of the n-th state is Eb(n) = Eb(1s) / (2*n - 1)^2.
# Calculate the denominator
denominator = (2 * n - 1)**2
# Calculate the final binding energy for the n=3 state
Eb_n = Eb_1s / denominator

print(f"Step 2: Calculate the binding energy for n = {n} using the 2D exciton model.")
print(f"The formula is: Eb(n) = Eb(1s) / (2*n - 1)^2")
print(f"For n = {n}:")
# Show the calculation with numbers substituted in
print(f"Eb({n}) = {Eb_1s} eV / (2 * {n} - 1)^2")
print(f"Eb({n}) = {Eb_1s} eV / ({2*n-1})^2")
print(f"Eb({n}) = {Eb_1s} eV / {denominator}")
print("-" * 30)
print(f"The Rydberg energy for n = {n} is: {Eb_n:.2f} eV")
