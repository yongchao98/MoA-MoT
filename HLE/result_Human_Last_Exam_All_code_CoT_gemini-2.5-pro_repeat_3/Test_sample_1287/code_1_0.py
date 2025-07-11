# Define monoisotopic masses of atoms and sodium ion
MASS_C = 12.0
MASS_H = 1.007825
MASS_N = 14.003074
MASS_O = 15.994915
MASS_NA = 22.989770

# 1. Start with the known m/z of the sodiated, permethylated, methyl-esterified A2G2S2 glycan.
mass_esterified_sodiated = 2901.4284

# 2. Calculate the neutral mass of this reference glycan by subtracting the mass of sodium.
mass_esterified_neutral = mass_esterified_sodiated - MASS_NA

# 3. Calculate the mass difference when a methyl ester group (-COOCH3) is replaced
#    by a dimethylated amide group (-CON(CH3)2).
# Mass of a methyl ester group: 2*C + 3*H + 2*O
mass_methyl_ester = 2 * MASS_C + 3 * MASS_H + 2 * MASS_O
# Mass of a dimethylated amide group: 3*C + 6*H + 1*N + 1*O
mass_dimethyl_amide = 3 * MASS_C + 6 * MASS_H + 1 * MASS_N + 1 * MASS_O

# The mass difference for one sialic acid modification
mass_diff_per_site = mass_dimethyl_amide - mass_methyl_ester

# 4. The glycan has two sialic acids, so the total mass difference is twice this value.
total_mass_diff = 2 * mass_diff_per_site

# 5. Calculate the neutral mass of the target amidated glycan.
mass_amidated_neutral = mass_esterified_neutral + total_mass_diff

# 6. Calculate the final m/z of the sodiated ion.
final_mz = mass_amidated_neutral + MASS_NA

# 7. Print the results.
# Since all three are isomers, they will have the same mass.
print("The glycans A2G(4)2S(3)2, A2G(4)S(3)S(6), and A2G(4)2S(6)2 are isomers and will have the same m/z.")
print("The expected m/z for the singly-sodiated ions is calculated as follows:\n")
print(f"1. Neutral mass of reference permethylated/esterified glycan:")
print(f"   {mass_esterified_sodiated:.4f} (M_ref+Na) - {MASS_NA:.4f} (Na) = {mass_esterified_neutral:.4f}\n")
print(f"2. Mass difference from converting two ester groups to dimethyl-amide groups:")
print(f"   2 * ({mass_dimethyl_amide:.4f} - {mass_methyl_ester:.4f}) = {total_mass_diff:.4f}\n")
print(f"3. Neutral mass of the target permethylated/amidated glycan:")
print(f"   {mass_esterified_neutral:.4f} + {total_mass_diff:.4f} = {mass_amidated_neutral:.4f}\n")
print(f"4. Final m/z of the sodiated ion [M+Na]+:")
print(f"   {mass_amidated_neutral:.4f} (M_final) + {MASS_NA:.4f} (Na) = {final_mz:.4f}\n")

print(f"Therefore, the expected m/z value for all three glycans is {final_mz:.4f}.")
<<<2927.4917>>>