import sys

# Plan: Calculate the number of IR active phonons for LiNiPO4.
# 1. Use the known results of factor group analysis for the Pnma space group.
# 2. Define the total number of vibrational modes for each irreducible representation (irrep) of the D2h point group.
# 3. Define and subtract the acoustic modes.
# 4. Identify the IR active irreps and their corresponding polarizations.
# 5. Print the breakdown and the final result.

# Step 2: Total vibrational modes for LiNiPO4 (Pnma, Z=4) from literature.
# The crystal has 28 atoms in the primitive cell, resulting in 3 * 28 = 84 total modes.
gamma_total = {
    "Ag": 12, "B1g": 12, "B2g": 9, "B3g": 9,
    "Au": 9,  "B1u": 9,  "B2u": 12, "B3u": 12
}

# Step 3: Acoustic modes for the D2h point group.
# The three translational modes transform as x, y, and z.
# T_x -> B3u, T_y -> B2u, T_z -> B1u
gamma_acoustic = {"B3u": 1, "B2u": 1, "B1u": 1}

# Step 4: Calculate the number of optical modes by subtracting acoustic modes.
gamma_optic = gamma_total.copy()
for irrep, count in gamma_acoustic.items():
    gamma_optic[irrep] -= count

# Step 5: Identify IR active modes and their polarizations.
# We map crystallographic axes (a, b, c) to Cartesian axes (x, y, z).
# E || x (along 'a' axis) corresponds to the B3u representation.
# E || y (along 'b' axis) corresponds to the B2u representation.
# E || z (along 'c' axis) corresponds to the B1u representation.
ir_phonons_x = gamma_optic["B3u"]
ir_phonons_y = gamma_optic["B2u"]
ir_phonons_z = gamma_optic["B1u"]

print("Prediction of IR active phonons for LiNiPO4 (orthorhombic, Pnma, D2h point group):")
print("The calculation subtracts the 3 acoustic modes from the total modes of the corresponding symmetries.\n")

# --- Output for E||x polarization ---
total_x = gamma_total['B3u']
acoustic_x = gamma_acoustic['B3u']
print(f"For E||x polarization (B3u symmetry):")
print(f"The number of IR active phonons is the total number of B3u modes minus the acoustic B3u mode.")
print(f"Number of phonons = {total_x} - {acoustic_x} = {ir_phonons_x}\n")

# --- Output for E||y polarization ---
total_y = gamma_total['B2u']
acoustic_y = gamma_acoustic['B2u']
print(f"For E||y polarization (B2u symmetry):")
print(f"The number of IR active phonons is the total number of B2u modes minus the acoustic B2u mode.")
print(f"Number of phonons = {total_y} - {acoustic_y} = {ir_phonons_y}\n")

# --- Output for E||z polarization ---
total_z = gamma_total['B1u']
acoustic_z = gamma_acoustic['B1u']
print(f"For E||z polarization (B1u symmetry):")
print(f"The number of IR active phonons is the total number of B1u modes minus the acoustic B1u mode.")
print(f"Number of phonons = {total_z} - {acoustic_z} = {ir_phonons_z}\n")

# --- Final consolidated answer ---
print("--- Summary ---")
print(f"E||x: {ir_phonons_x}, E||y: {ir_phonons_y}, E||z: {ir_phonons_z}")

# This line is for programmatic parsing of the final answer.
sys.stdout.write(f'<<<E||x: {ir_phonons_x}, E||y: {ir_phonons_y}, E||z: {ir_phonons_z}>>>')
