# --- 1. Define constants and initial composition ---

# Monoisotopic atomic masses
MASS_C = 12.000000
MASS_H = 1.007825
MASS_O = 15.994915
MASS_N = 14.003074
MASS_Na_ion = 22.989221  # Mass of Na+ ion ([23Na] - e-)

# Glycan composition (standard biantennary, disialylated N-glycan: A2G2S2)
composition = {
    'Man': {'formula': {'C': 6, 'H': 12, 'O': 6, 'N': 0}, 'count': 3},
    'GlcNAc': {'formula': {'C': 8, 'H': 15, 'O': 6, 'N': 1}, 'count': 4},
    'Gal': {'formula': {'C': 6, 'H': 12, 'O': 6, 'N': 0}, 'count': 2},
    'Neu5Ac': {'formula': {'C': 11, 'H': 19, 'O': 9, 'N': 1}, 'count': 2}
}

glycan_names = ["A2G(4)2S(3)2", "A2G(4)S(3)S(6)", "A2G(4)2S(6)2"]

# --- 2. Calculate the formula of the base glycan polymer ---

base_formula = {'C': 0, 'H': 0, 'O': 0, 'N': 0}
total_monomers = 0

# Sum the formulas of all monomers
for sugar in composition.values():
    total_monomers += sugar['count']
    for atom, count in sugar['formula'].items():
        base_formula[atom] += count * sugar['count']

# Subtract water molecules (H2O) for each glycosidic bond formed
# For a polymer of n monomers, there are n-1 bonds. The glycan is not cyclic.
num_bonds = total_monomers - 1
base_formula['H'] -= 2 * num_bonds
base_formula['O'] -= 1 * num_bonds

# --- 3. Apply chemical modifications to the formula ---

# Step A: Amidation of 2 sialic acids (COOH -> CONH2)
# Change per site: -O, +N, +H
num_amide_reactions = 2
modified_formula = base_formula.copy()
modified_formula['O'] -= 1 * num_amide_reactions
modified_formula['N'] += 1 * num_amide_reactions
modified_formula['H'] += 1 * num_amide_reactions

# Step B: Permethylation (replace H with CH3 on all free -OH groups)
# Net change per site: +CH2
# For a standard A2G2S2 N-glycan, there are 29 methylation sites.
num_methylation_sites = 29
modified_formula['C'] += 1 * num_methylation_sites
modified_formula['H'] += 2 * num_methylation_sites

# --- 4. Calculate the neutral mass and final m/z ---

# Calculate the monoisotopic mass of the final neutral molecule
neutral_mass = (modified_formula['C'] * MASS_C +
                modified_formula['H'] * MASS_H +
                modified_formula['O'] * MASS_O +
                modified_formula['N'] * MASS_N)

# Calculate the m/z for the singly-sodiated ion [M+Na]+
mz_sodiated = neutral_mass + MASS_Na_ion

# --- 5. Print the results ---

print("Since all three glycans are isomers that undergo the same chemical reactions,")
print("they will all have the same final mass-to-charge ratio.\n")
print("Calculation Steps:")
print(f"1. Base glycan polymer formula: C{base_formula['C']}H{base_formula['H']}N{base_formula['N']}O{base_formula['O']}")
print(f"2. After amidation & permethylation, final formula: C{modified_formula['C']}H{modified_formula['H']}N{modified_formula['N']}O{modified_formula['O']}")
print("\nCalculating m/z for the [M+Na]+ ion:")
print(f"Final Neutral Mass (M) = {neutral_mass:.4f} Da")
print(f"Sodium Ion Mass (Na+) = {MASS_Na_ion:.4f} Da")
print(f"Final Equation: {neutral_mass:.4f} (M) + {MASS_Na_ion:.4f} (Na+) = {mz_sodiated:.4f} m/z\n")

print("Expected singly-sodiated mass [M+Na]+ for the analyzed glycans:")
for name in glycan_names:
    print(f"- {name}: {mz_sodiated:.4f} m/z")

# Final answer in specified format
final_answer_val = f"{mz_sodiated:.4f}"
print(f"\n<<<Solution_End>>>")
print(f"The calculated m/z value is {final_answer_val}")