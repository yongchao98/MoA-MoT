# Atomic masses (integer masses for this level of analysis)
mass_C = 12
mass_H = 1
mass_Cl = 35

# Proposed Molecular Formula: C12H11Cl3
molecular_ion_mass = 12 * mass_C + 11 * mass_H + 3 * mass_Cl
print(f"Proposed Molecular Ion [C12H11Cl3]+ calculation:")
print(f"12 * {mass_C} (C) + 11 * {mass_H} (H) + 3 * {mass_Cl} (Cl) = {molecular_ion_mass}")
print(f"This matches the observed M+ peak at m/z 260.\n")

# Base Peak: Loss of one Cl
base_peak_mass = molecular_ion_mass - mass_Cl
print(f"Base Peak [M-Cl]+ calculation:")
print(f"{molecular_ion_mass} - {mass_Cl} (Cl) = {base_peak_mass}")
print(f"This matches the observed base peak at m/z 225.\n")

# Second Fragmentation: Loss of second Cl
second_fragment_mass = base_peak_mass - mass_Cl
print(f"Second Fragment [M-2Cl]+ calculation:")
print(f"{base_peak_mass} - {mass_Cl} (Cl) = {second_fragment_mass}")
print(f"This matches the observed peak at m/z 190.\n")

# Third Fragmentation: Loss of third Cl
third_fragment_mass = second_fragment_mass - mass_Cl
print(f"Third Fragment [M-3Cl]+ calculation:")
print(f"{second_fragment_mass} - {mass_Cl} (Cl) = {third_fragment_mass}")
print(f"This is the [C12H11]+ fragment and matches the observed peak at m/z 155.\n")

# Key Structural Fragments
# Trichlorocyclopropenyl cation
c3cl3_fragment_mass = 3 * mass_C + 3 * mass_Cl
print(f"Key Fragment [C3Cl3]+ calculation:")
print(f"3 * {mass_C} (C) + 3 * {mass_Cl} (Cl) = {c3cl3_fragment_mass}")
print(f"This matches the observed peak with a 3-Cl isotope pattern at m/z 141.\n")

# Mesityl cation
c9h11_fragment_mass = 9 * mass_C + 11 * mass_H
print(f"Key Fragment [C9H11]+ (Mesityl) calculation:")
print(f"9 * {mass_C} (C) + 11 * {mass_H} (H) = {c9h11_fragment_mass}")
print(f"This matches the observed peak at m/z 119.\n")

# Final Conclusion
compound_name = "1,2,3-trichloro-3-(2,4,6-trimethylphenyl)cycloprop-1-ene"
print(f"The evidence strongly suggests the compound is:")
print(compound_name)