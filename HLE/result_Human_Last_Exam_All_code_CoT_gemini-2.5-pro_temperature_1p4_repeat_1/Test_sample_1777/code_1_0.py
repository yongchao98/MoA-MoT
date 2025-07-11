import sys

# --- Given Parameters ---
Eg = 3.0  # Band gap in eV
E_peak_1s = 1.0  # 1s exciton resonance peak in eV
n = 3     # Target principal quantum number

# --- Step 1: Calculate the 1s Exciton Binding Energy (E_b_1s) ---
print("Step 1: Calculate the 1s exciton binding energy.")
print(f"E_binding_1s = Band Gap - E_peak_1s")
print(f"E_binding_1s = {Eg} eV - {E_peak_1s} eV")
E_b_1s = Eg - E_peak_1s
print(f"E_binding_1s = {E_b_1s} eV")
print("-" * 20)

# --- Step 2: Calculate the Effective Rydberg Energy (Ry*) for a 2D system ---
# The formula for 2D exciton binding energy is: E_binding_n = Ry* / (n - 0.5)^2
# For n=1, this is E_binding_1s = Ry* / (1 - 0.5)^2 = Ry* / 0.25 = 4 * Ry*
# Therefore, Ry* = E_binding_1s / 4
print("Step 2: Calculate the effective Rydberg energy (Ry*).")
print(f"For a 2D system, Ry* = E_binding_1s / (1 - 0.5)^2 = E_binding_1s / 4")
Ry_star = E_b_1s / 4
print(f"Ry* = {E_b_1s} eV / 4")
print(f"Ry* = {Ry_star} eV")
print("-" * 20)

# --- Step 3: Calculate the Binding Energy for n=3 (E_b_3) ---
print(f"Step 3: Calculate the Rydberg (binding) energy for n={n}.")
print(f"E_binding_{n} = Ry* / ({n} - 0.5)^2")
denominator = (n - 0.5) ** 2
print(f"E_binding_{n} = {Ry_star} eV / ({n - 0.5})^2")
print(f"E_binding_{n} = {Ry_star} eV / {denominator}")
E_b_3 = Ry_star / denominator
# Use sys.stdout.write to prevent adding an extra newline, making the final answer tag appear on its own line.
sys.stdout.write(f"\nThe Rydberg energy for n = {n} is {E_b_3:.2f} eV.\n")

sys.stdout.write(f"<<<{E_b_3:.2f}>>>")