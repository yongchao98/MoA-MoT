# Given parameters
Eg = 3.0  # Band gap in eV
E_1s = 1.0 # Resonance peak for 1s exciton in eV
n1 = 1      # Principal quantum number for the 1s state
n3 = 3      # Target principal quantum number

print("Step 1: Calculate the binding energy of the 1s exciton (Eb_1s).")
# The exciton energy is the band gap minus its binding energy: E_n = Eg - Eb_n
# Therefore, the binding energy is: Eb_n = Eg - E_n
Eb_1s = Eg - E_1s
print(f"The binding energy for n={n1} is Eb({n1}) = Eg - E_({n1}s) = {Eg} eV - {E_1s} eV = {Eb_1s} eV\n")


print("Step 2: Calculate the effective Rydberg energy (R_y*) for the 2D system.")
# For a 2D system, the binding energy is given by: Eb_n = R_y* / (n - 0.5)^2
# We can find R_y* from the 1s state: R_y* = Eb_1s * (n - 0.5)^2
R_y_star = Eb_1s * (n1 - 0.5)**2
print(f"The effective Rydberg energy R_y* = Eb({n1}) * ({n1} - 0.5)^2 = {Eb_1s} eV * {((n1 - 0.5)**2)} = {R_y_star} eV\n")


print("Step 3: Calculate the binding energy for the n=3 state (Eb_3).")
# This is the value requested by the user, often called the 'Rydberg energy for n=3'.
Eb_3 = R_y_star / (n3 - 0.5)**2
print(f"The binding energy for n={n3}, Eb({n3}), is calculated as:")
print(f"Eb({n3}) = R_y* / ({n3} - 0.5)^2 = {R_y_star} eV / ({(n3 - 0.5)**2}) = {Eb_3:.3f} eV")

print("\nThe Rydberg energy (binding energy) for n = 3 is the final result.")
print(f"Final Answer: {Eb_3:.3f} eV")