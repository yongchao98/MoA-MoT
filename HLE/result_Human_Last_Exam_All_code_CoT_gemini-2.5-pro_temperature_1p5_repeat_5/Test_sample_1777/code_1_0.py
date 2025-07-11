# Define the given parameters
Eg = 3.0  # Band gap in eV
E_ex_1s = 1.0  # 1s exciton resonance peak in eV
n = 3  # Target quantum number

# Step 1: Calculate the binding energy of the 1s exciton state (Eb_1s)
Eb_1s = Eg - E_ex_1s
print(f"The binding energy of the 1s state is calculated by subtracting the resonance peak energy from the band gap:")
print(f"Eb(1s) = {Eg} eV - {E_ex_1s} eV = {Eb_1s} eV\n")

# Step 2: Calculate the exciton Rydberg energy (Ry) for the 2D system
# For a 2D system, Eb(n) = Ry / (n - 0.5)^2
# So, for n=1, Ry = Eb(1s) * (1 - 0.5)^2
n_1s = 1
Ry = Eb_1s * (n_1s - 0.5)**2
print(f"The exciton Rydberg energy (Ry) is found using the 2D exciton formula for n=1:")
print(f"Ry = Eb(1s) * ({n_1s} - 0.5)^2 = {Eb_1s} eV * {((n_1s - 0.5)**2)} = {Ry} eV\n")

# Step 3: Calculate the Rydberg energy (binding energy) for the n=3 state
Eb_ns = Ry / (n - 0.5)**2
print(f"The Rydberg energy for n={n} is the binding energy of that state, calculated as:")
print(f"Eb(n={n}) = Ry / ({n} - 0.5)^2")
print(f"Eb(n={n}) = {Ry} eV / ({n - 0.5})^2 = {Eb_ns:.3f} eV")

# Final answer
final_answer = round(Eb_ns, 3)