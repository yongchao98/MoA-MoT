import sys

# Define given values
Eg = 3.0  # Band gap in eV
E_1s = 1.0 # 1s exciton resonance peak in eV
n = 3      # Principal quantum number for the target state

# --- Step 1: Calculate the binding energy of the 1s exciton ---
# The energy of the 1s exciton state is the band gap minus its binding energy.
# E_1s = Eg - E_b,1s  =>  E_b,1s = Eg - E_1s
Eb_1s = Eg - E_1s
print("Step 1: Calculate the ground state (1s) exciton binding energy (Eb_1s).")
print(f"Eb_1s = Eg - E_1s = {Eg} eV - {E_1s} eV = {Eb_1s} eV\n")

# --- Step 2: Calculate the exciton Rydberg constant (Rx) for a 2D system ---
# For a 2D system, the binding energy series is: Eb,n = Rx / (n - 1/2)^2
# For the ground state (n=1), this gives: Eb,1s = Rx / (1 - 0.5)^2 = 4 * Rx
# Therefore, Rx = Eb,1s / 4
Rx = Eb_1s / 4.0
print("Step 2: Calculate the exciton Rydberg constant (Rx) for a 2D system.")
print("Eb,1s = Rx / (1 - 0.5)^2 = 4 * Rx")
print(f"Rx = Eb_1s / 4 = {Eb_1s} eV / 4 = {Rx} eV\n")

# --- Step 3: Calculate the Rydberg energy for n=3 ---
# We use the 2D formula again for n=3.
Eb_n = Rx / (n - 0.5)**2
print(f"Step 3: Calculate the Rydberg energy (binding energy) for n={n}.")
# Output the final equation with all numbers as requested
print(f"The final equation is: Eb,{n} = Rx / ({n} - 0.5)^2")
print(f"Plugging in the numbers: Eb,{n} = {Rx} / ({n} - 0.5)^2 = {Eb_n:.3f} eV")

# Format and write the final answer to stdout for grading
# The sys.stdout.write is used to ensure no extra characters are printed
final_answer = f"<<<{Eb_n:.3f}>>>"
sys.stdout.write(final_answer)
