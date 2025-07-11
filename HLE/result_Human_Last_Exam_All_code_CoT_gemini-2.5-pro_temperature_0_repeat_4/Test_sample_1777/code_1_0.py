import math

# Given parameters
Eg = 3.0  # Band gap in eV
E_1s = 1.0 # 1s exciton resonance peak in eV
n_target = 3 # Target principal quantum number

# Step 1: Calculate the binding energy of the 1s exciton state
# E_b(n) = Eg - E_n
Eb_1s = Eg - E_1s
print(f"Step 1: Calculate the binding energy of the 1s state.")
print(f"E_b(1s) = E_g - E_1s = {Eg} eV - {E_1s} eV = {Eb_1s} eV\n")

# Step 2: Determine the effective Rydberg energy (Ry) for the 2D system
# For a 2D exciton, E_b(n) = Ry / (n - 0.5)^2
# So, for n=1, E_b(1s) = Ry / (1 - 0.5)^2 = Ry / 0.25 = 4 * Ry
# Therefore, Ry = E_b(1s) / 4
Ry = Eb_1s / 4.0
print(f"Step 2: Determine the effective Rydberg energy (Ry).")
print(f"From E_b(1s) = Ry / (1 - 0.5)^2, we get Ry = E_b(1s) / 4 = {Eb_1s} eV / 4 = {Ry} eV\n")

# Step 3: Calculate the binding energy for the n=3 state
# E_b(3) = Ry / (3 - 0.5)^2
denominator = (n_target - 0.5)**2
Eb_target = Ry / denominator
print(f"Step 3: Calculate the binding energy for n = {n_target}.")
print("The final equation is:")
print(f"E_b({n_target}) = Ry / ({n_target} - 0.5)^2 = {Ry} eV / ({n_target - 0.5})^2 = {Eb_target:.3f} eV")
