# --- Atomic and Molecular Masses (monoisotopic) ---
mass_H = 1.007825
mass_C = 12.000000
mass_N = 14.003074
mass_O = 15.994915
mass_Na = 22.989770
mass_H2O = 2 * mass_H + mass_O

# --- Glycan Residue Composition ---
# Hexose (e.g., Galactose, Mannose): C6H10O5
mass_hex_residue = 6 * mass_C + 10 * mass_H + 5 * mass_O
# N-Acetylhexosamine (e.g., GlcNAc): C8H13NO5
mass_hexnac_residue = 8 * mass_C + 13 * mass_H + 1 * mass_N + 5 * mass_O
# N-Acetylneuraminic acid (Sialic acid, NeuAc): C11H17NO8
mass_neuac_residue = 11 * mass_C + 17 * mass_H + 1 * mass_N + 8 * mass_O

# --- Glycan Composition ---
# A2G2S2: 5 Hexose, 4 HexNAc, 2 NeuAc
num_hex = 5
num_hexnac = 4
num_neuac = 2

# 1. Calculate the initial mass of the native glycan (A2G2S2)
# The mass is the sum of residues plus one water molecule for the reducing end.
mass_initial_glycan = (num_hex * mass_hex_residue +
                       num_hexnac * mass_hexnac_residue +
                       num_neuac * mass_neuac_residue +
                       mass_H2O)

# 2. Calculate mass after amidation
# Reaction: R-COOH + NH3 -> R-CONH2 + H2O. Net change is -OH + NH2 per reaction.
# There are two sialic acids, so two reactions.
mass_change_amidation = -mass_O - mass_H + mass_N + 2 * mass_H # -OH + NH2
mass_amidated_glycan = mass_initial_glycan + 2 * mass_change_amidation

# 3. Calculate mass after permethylation
# We need to count the number of replaceable H atoms on O and N.
# On the amidated glycan, this corresponds to:
# - OH groups: 31
# - N-H groups of N-acetyl functions: 6
# - N-H groups of the two new primary amides (-CONH2): 2 * 2 = 4
num_hydroxyls = 31
num_n_acetyl = 6
num_amide_hydrogens = 4
total_methylations = num_hydroxyls + num_n_acetyl + num_amide_hydrogens

# Mass change for permethylation: replace H with CH3. Net change is addition of CH2.
mass_change_per_methylation = mass_C + 2 * mass_H
total_mass_change_methylation = total_methylations * mass_change_per_methylation
mass_derivatized_neutral = mass_amidated_glycan + total_mass_change_methylation

# 4. Calculate the final m/z for the sodiated ion [M+Na]+
final_mz = mass_derivatized_neutral + mass_Na

# --- Print the results step-by-step ---
print("Calculation of the expected m/z for the derivatized glycans.")
print("-" * 60)
print(f"Initial A2G2S2 glycan mass: {mass_initial_glycan:.4f} Da")
print(f"Mass after amidation of 2 sialic acids: {mass_amidated_glycan:.4f} Da")
print(f"Number of methylation sites (31 OH + 6 NAc + 2 CONH2): {total_methylations}")
print(f"Mass of the fully derivatized neutral glycan: {mass_derivatized_neutral:.4f} Da")
print("-" * 60)
print("Final m/z Calculation for [M+Na]+ ion:")
print(f"Mass of derivatized glycan + Mass of Sodium = Final m/z")
print(f"{mass_derivatized_neutral:.4f} Da + {mass_Na:.4f} Da = {final_mz:.4f} Da")
print("\nSince the atomic composition is identical regardless of the sialic acid linkage,")
print("the expected m/z value is the same for all three glycans.")
print(f"Mass for A2G(4)2S(3)2: {final_mz:.4f}")
print(f"Mass for A2G(4)S(3)S(6): {final_mz:.4f}")
print(f"Mass for A2G(4)2S(6)2: {final_mz:.4f}")
print(f"<<<{final_mz:.4f}>>>")