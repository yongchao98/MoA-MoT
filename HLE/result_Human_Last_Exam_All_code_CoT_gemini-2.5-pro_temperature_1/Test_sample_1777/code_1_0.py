import sys

# Given parameters
E_g = 3.0  # Band gap in eV
E_1s = 1.0 # 1s exciton resonance peak in eV
n = 3      # Principal quantum number for the target state

# --- Step 1: Calculate the 1s exciton binding energy (E_b_1s) ---
# The binding energy is the difference between the band gap and the exciton state energy.
E_b_1s = E_g - E_1s
print(f"Step 1: Calculate the binding energy of the 1s exciton (E_b_1s).")
print(f"Equation: E_b_1s = E_g - E_1s")
print(f"Calculation: E_b_1s = {E_g} eV - {E_1s} eV = {E_b_1s} eV\n")

# --- Step 2: Calculate the effective Rydberg energy (Ry_star) for a 2D system ---
# For a 2D exciton, the binding energy is E_b,n = Ry* / (n - 0.5)^2.
# We can find Ry* from the 1s state: Ry* = E_b_1s * (1 - 0.5)^2.
n1 = 1
n1_term = (n1 - 0.5)**2
Ry_star = E_b_1s * n1_term
print(f"Step 2: Calculate the effective Rydberg energy (Ry*).")
print(f"Equation: Ry* = E_b_1s * (n - 0.5)^2 for n=1")
print(f"Calculation: Ry* = {E_b_1s} eV * ({n1} - 0.5)^2 = {Ry_star} eV\n")

# --- Step 3: Calculate the Rydberg energy (binding energy) for n=3 ---
# Use the calculated Ry_star to find the binding energy for the n=3 state.
n3_term = (n - 0.5)**2
E_b_3 = Ry_star / n3_term
print(f"Step 3: Calculate the Rydberg energy for n = {n}.")
print(f"Equation: E_b_{n} = Ry* / (n - 0.5)^2")
print(f"Calculation: E_b_{n} = {Ry_star} eV / ({n} - 0.5)^2 = {E_b_3:.2f} eV\n")

# --- Final Answer ---
final_answer = E_b_3
# Use sys.stdout to ensure the final answer is the last line printed
sys.stdout.write(f"<<<{final_answer:.2f}>>>\n")