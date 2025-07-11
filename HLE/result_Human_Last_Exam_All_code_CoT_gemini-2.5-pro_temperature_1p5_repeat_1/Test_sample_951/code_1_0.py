# Step 1 & 2: Define the system and its properties
# Crystal: LiNiPO4 (olivine structure)
# Space Group: Pnma (No. 62)
# Point Group: D2h

# Step 3: Decomposition of total vibrational modes into irreducible representations
# This result is derived from group theory analysis, often using specialized software
# or tables based on the crystal's symmetry and atomic positions.
# The irreps for D2h are Ag, B1g, B2g, B3g, Au, B1u, B2u, B3u.
total_modes = {
    "Ag": 12,
    "B1g": 9,
    "B2g": 12,
    "B3g": 9,
    "Au": 9,
    "B1u": 12,
    "B2u": 9,
    "B3u": 12,
}

# Step 4: Identify IR-active modes and subtract acoustic modes
# IR-active modes are of ungerade ('u') symmetry: B1u, B2u, B3u.
# Acoustic modes also have 'u' symmetry: one each for B1u, B2u, B3u.
# To get the number of IR-active OPTICAL phonons, we subtract 1 from each.
ir_optical_modes = {
    "B1u": total_modes["B1u"] - 1,
    "B2u": total_modes["B2u"] - 1,
    "B3u": total_modes["B3u"] - 1,
}

# Step 5: Assign polarizations and calculate final numbers
# For the D2h point group in a standard axis setting:
# E||x excites B3u modes.
# E||y excites B2u modes.
# E||z excites B1u modes.
num_phonons_Ex = ir_optical_modes["B3u"]
num_phonons_Ey = ir_optical_modes["B2u"]
num_phonons_Ez = ir_optical_modes["B1u"]

# Final Output: Print the result in the requested format
# The output will show each number in the final equation.
print(f"Total B1u modes: {total_modes['B1u']}")
print(f"Total B2u modes: {total_modes['B2u']}")
print(f"Total B3u modes: {total_modes['B3u']}")
print("\nSubtracting 1 acoustic mode from each IR-active symmetry:")
print(f"Optical B1u (E||z) = {total_modes['B1u']} - 1 = {num_phonons_Ez}")
print(f"Optical B2u (E||y) = {total_modes['B2u']} - 1 = {num_phonons_Ey}")
print(f"Optical B3u (E||x) = {total_modes['B3u']} - 1 = {num_phonons_Ex}")
print("\nFinal prediction:")
print(f"E||x: {num_phonons_Ex}, E||y: {num_phonons_Ey}, E||z: {num_phonons_Ez}")
