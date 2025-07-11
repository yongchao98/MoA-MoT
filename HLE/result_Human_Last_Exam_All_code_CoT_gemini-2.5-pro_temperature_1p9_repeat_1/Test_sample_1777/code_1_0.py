# Given parameters
Eg = 3.0  # Band gap in eV
E_peak_1s = 1.0  # Resonance peak for 1s exciton in eV
n = 3 # The target principal quantum number

# --- Step 1: Calculate the 1s Exciton Binding Energy ---
# The binding energy is the difference between the band gap and the resonance peak energy.
Eb_1s = Eg - E_peak_1s
print(f"Step 1: Calculate the 1s exciton binding energy (Eb_1s).")
print(f"Eb_1s = Band Gap - 1s Resonance Peak = {Eg} eV - {E_peak_1s} eV = {Eb_1s} eV\n")

# --- Step 2: Use the 2D exciton model to find the n=3 binding energy ---
# The binding energy for a 2D exciton is E_b(n) = R_y' / (n - 1/2)^2.
# We can find the ratio E_b(3)/E_b(1).
# E_b(3) / E_b(1) = [(1 - 0.5)^2] / [(3 - 0.5)^2] = 0.5^2 / 2.5^2 = 1 / 25
# So, E_b(3) = E_b(1) / 25
denominator = ((n - 0.5)**2) / ((1 - 0.5)**2)
Eb_3 = Eb_1s / denominator

# --- Step 3: Print the final calculation and result ---
print(f"Step 2: Calculate the n=3 binding energy (Eb_3).")
print(f"For a 2D material, the binding energy Eb_n scales as 1/(n-0.5)^2.")
print(f"Therefore, Eb_3 = Eb_1s / ( (3-0.5)^2 / (1-0.5)^2 )")
print(f"Eb_3 = {Eb_1s} eV / ( (2.5)^2 / (0.5)^2 )")
print(f"Eb_3 = {Eb_1s} eV / ( {2.5**2} / {0.5**2} )")
print(f"Eb_3 = {Eb_1s} eV / {denominator}")
print(f"\nThe final Rydberg energy (binding energy) for n = 3 is: {Eb_3:.2f} eV")

# The final numerical answer
final_answer = Eb_3