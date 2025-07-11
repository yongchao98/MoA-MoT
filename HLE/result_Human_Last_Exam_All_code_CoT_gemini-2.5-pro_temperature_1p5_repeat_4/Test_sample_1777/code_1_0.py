# Given parameters
Eg = 3.0  # Band gap in eV
E1s = 1.0 # 1s exciton resonance peak in eV
n = 3     # Target principal quantum number

# --- Step 1: Calculate the binding energy of the 1s state (ground state) ---
# The exciton energy is the band gap minus its binding energy: E_n = Eg - EB_n
# So, the binding energy is: EB_n = Eg - E_n
EB_1s = Eg - E1s

# --- Step 2: Determine the 2D effective Rydberg energy (Ry_2D) ---
# For a 2D exciton, the binding energy series is EB_n = Ry_2D / (n - 0.5)^2
# For n=1, EB_1s = Ry_2D / (1 - 0.5)^2. We solve for Ry_2D.
Ry_2D = EB_1s * (1 - 0.5)**2

# --- Step 3: Calculate the binding energy for n=3 ---
# The "Rydberg energy for n=3" refers to the binding energy of the n=3 state.
EB_n = Ry_2D / (n - 0.5)**2

# --- Final Output ---
print(f"The binding energy for the n={n} state is calculated using the 2D exciton model:")
print(f"EB_n = Ry_2D* / (n - 1/2)^2")
print("\nFirst, we find the 2D effective Rydberg energy (Ry_2D*):")
print(f"Ry_2D* = (Eg - E1s) * (1 - 0.5)^2 = ({Eg} - {E1s}) * {(1 - 0.5)**2} = {Ry_2D} eV")
print("\nNow, we calculate the binding energy for n=3:")
# Final equation with numbers
equation = f"EB_3 = {Ry_2D:.2f} eV / ({n} - 0.5)^2 = {EB_n:.2f} eV"
print(equation)