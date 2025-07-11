#
# Plan:
# 1. Calculate the binding energy of the 1s exciton (E_b_1s) using the band gap (E_g) and the 1s resonance peak (E_1s_peak).
#    E_b_1s = E_g - E_1s_peak
# 2. Determine the effective Rydberg constant (R_y) using the formula for 2D excitons: E_b_n = R_y / (n - 0.5)^2.
#    For n=1, this is E_b_1s = R_y / (1 - 0.5)^2, which gives R_y = E_b_1s * (0.5)^2.
#    Hold on, the formula is E_b,n = R_y / (n-0.5)^2.
#    So, E_b,1s = R_y / (1-0.5)^2 = R_y / (0.5)^2 = R_y / 0.25 = 4*R_y.
#    Therefore, R_y = E_b_1s / 4.
# 3. Calculate the binding energy for the n=3 state (E_b_3), which is the quantity requested.
#    E_b_3 = R_y / (3 - 0.5)^2
# 4. Print the final result, showing the numbers used in the final calculation.
#

# Given values in eV
band_gap_Eg = 3.0
resonance_peak_E1 = 1.0

# State we are interested in
n = 3

# Step 1: Calculate the binding energy of the 1s exciton
binding_energy_E_b1 = band_gap_Eg - resonance_peak_E1

# Step 2: Calculate the effective Rydberg constant (R_y) for the 2D system
# E_b,1s = R_y / (1 - 0.5)^2 = 4 * R_y
rydberg_constant_Ry = binding_energy_E_b1 / 4.0

# Step 3: Calculate the binding energy for the n=3 state
n_minus_half = n - 0.5
binding_energy_E_b3 = rydberg_constant_Ry / (n_minus_half)**2

# Step 4: Print the results including the final equation
print(f"Given Band Gap (E_g): {band_gap_Eg} eV")
print(f"Given 1s Exciton Resonance Peak (E_1s): {resonance_peak_E1} eV")
print("-" * 30)
print(f"Step 1: Calculate 1s Binding Energy (E_b,1s)")
print(f"E_b,1s = E_g - E_1s = {band_gap_Eg} eV - {resonance_peak_E1} eV = {binding_energy_E_b1} eV")
print("-" * 30)
print(f"Step 2: Calculate Effective Rydberg Constant (R_y)")
print(f"From E_b,1s = R_y / (1 - 0.5)^2  =>  R_y = E_b,1s / 4")
print(f"R_y = {binding_energy_E_b1} eV / 4 = {rydberg_constant_Ry} eV")
print("-" * 30)
print(f"Step 3: Calculate the Rydberg energy for n={n}")
print(f"The formula is E_b,n = R_y / (n - 0.5)^2")
print(f"E_b,{n} = {rydberg_constant_Ry} eV / ({n} - 0.5)^2")
print(f"E_b,{n} = {rydberg_constant_Ry} eV / ({n_minus_half})^2")
print(f"E_b,{n} = {rydberg_constant_Ry} eV / {n_minus_half**2}")
print(f"Final Answer: {binding_energy_E_b3:.2f} eV")
