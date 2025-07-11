import math

# Given parameters
Eg = 3.0  # Band gap in eV
E_1s = 1.0  # 1s exciton resonance peak in eV
n1 = 1  # Principal quantum number for the 1s state
n3 = 3  # Principal quantum number for the target state

# Step 1: Calculate the effective Rydberg energy (Ry) from the 1s exciton data.
# The energy of an exciton state in 2D is E_n = Eg - Ry / (n - 0.5)^2
# So, Ry = (Eg - E_n) * (n - 0.5)^2
Ry = (Eg - E_1s) * (n1 - 0.5)**2

# Step 2: Calculate the binding energy for the n=3 state.
# The question "Rydberg energy for n=3" refers to the binding energy of this state.
# E_binding(n) = Ry / (n - 0.5)^2
E_binding_3 = Ry / (n3 - 0.5)**2

# Step 3: Print the results and the equations with the numbers plugged in.
print("Step 1: Calculate the effective Rydberg energy (Ry) for the material.")
print(f"The formula for the 1s exciton energy is: E_1s = Eg - Ry / (n1 - 0.5)^2")
print(f"Rearranging to solve for Ry: Ry = (Eg - E_1s) * (n1 - 0.5)^2")
print(f"Plugging in the values: Ry = ({Eg} - {E_1s}) * ({n1} - 0.5)^2 = {Ry} eV")
print("\nStep 2: Calculate the binding energy for the n=3 state.")
print(f"The formula for the binding energy is: E_binding(n) = Ry / (n - 0.5)^2")
print(f"Plugging in the values for n = {n3}:")
print(f"E_binding({n3}) = {Ry} / ({n3} - 0.5)^2 = {E_binding_3:.2f} eV")

print(f"\nThe Rydberg energy for n = 3 is {E_binding_3:.2f} eV.")