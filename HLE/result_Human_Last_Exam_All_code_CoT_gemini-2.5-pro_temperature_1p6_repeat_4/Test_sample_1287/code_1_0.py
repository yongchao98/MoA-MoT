import collections

def calculate_mass(formula: dict, atomic_masses: dict) -> float:
    """Calculates the monoisotopic mass from a formula dictionary."""
    mass = 0.0
    for element, count in formula.items():
        mass += count * atomic_masses[element]
    return mass

# Monoisotopic masses of relevant elements and sodium
ATOMIC_MASSES = {
    'C': 12.0000000,
    'H': 1.0078250,
    'O': 15.9949146,
    'N': 14.0030740,
    'Na': 22.9897700,
}

# --- Step 1: Determine the formula of the native glycan ---
# Composition: Man3-GlcNAc4-Gal2-NeuAc2 (Total 11 residues)
# We calculate this from the sum of the full monosaccharides minus water for each of the 10 glycosidic bonds.
monosaccharide_formulas = {
    'Man': {'C': 6, 'H': 12, 'O': 6},
    'Gal': {'C': 6, 'H': 12, 'O': 6},
    'GlcNAc': {'C': 8, 'H': 15, 'N': 1, 'O': 6},
    'NeuAc': {'C': 11, 'H': 19, 'N': 1, 'O': 9},
    'H2O': {'H': 2, 'O': 1}
}
glycan_composition_counts = {'Man': 3, 'GlcNAc': 4, 'Gal': 2, 'NeuAc': 2}

native_glycan_formula = collections.defaultdict(int)
for name, count in glycan_composition_counts.items():
    for element, atom_count in monosaccharide_formulas[name].items():
        native_glycan_formula[element] += count * atom_count

# Subtract 10 water molecules for the 10 glycosidic bonds
num_bonds = sum(glycan_composition_counts.values()) - 1
for element, atom_count in monosaccharide_formulas['H2O'].items():
    native_glycan_formula[element] -= num_bonds * atom_count

# --- Step 2: Determine the formula of the final derivatized glycan ---
# The two reactions are:
# 1. Amidation: 2x sialic acid (-COOH -> -CONH2). Per site change: -O, +N, +H
# 2. Permethylation: All -OH and -NH groups are methylated (-H -> -CH3). Per site change: +C, +H2

# Calculate amidation changes
amidated_glycan_formula = native_glycan_formula.copy()
num_amidations = 2
amidated_glycan_formula['O'] -= num_amidations * 1
amidated_glycan_formula['N'] += num_amidations * 1
amidated_glycan_formula['H'] += num_amidations * 1

# Calculate number of methylation sites on the amidated glycan
# Sites = (sum of available sites on free monosaccharides) - 2 * (number of bonds)
# Free monosaccharide sites: Man(5), Gal(5), GlcNAc(5), Amidated-NeuAc(8)
sites_on_free_monosaccharides = (
    glycan_composition_counts['Man'] * 5 +
    glycan_composition_counts['Gal'] * 5 +
    glycan_composition_counts['GlcNAc'] * 5 +
    glycan_composition_counts['NeuAc'] * 8  # 5 OH + 1 original NH + 2 new amide NH
)
num_methyl_sites = sites_on_free_monosaccharides - 2 * num_bonds

# Calculate permethylation changes
final_glycan_formula = amidated_glycan_formula.copy()
final_glycan_formula['C'] += num_methyl_sites * 1
final_glycan_formula['H'] += num_methyl_sites * 2

# --- Step 3: Calculate the m/z of the [M+Na]+ ion ---
mass_M = calculate_mass(final_glycan_formula, ATOMIC_MASSES)
mass_Na = ATOMIC_MASSES['Na']
mass_M_plus_Na = mass_M + mass_Na

# --- Step 4: Print the results ---
print("All three glycans (A2G(4)2S(3)2, A2G(4)S(3)S(6), and A2G(4)2S(6)2) have the same composition and undergo identical reactions, resulting in the same final mass.")
print("\nCalculation Steps:")
print(f"1. Native Glycan Formula: C{native_glycan_formula['C']}H{native_glycan_formula['H']}N{native_glycan_formula['N']}O{native_glycan_formula['O']}")
print(f"2. After Amidation (x2): C{amidated_glycan_formula['C']}H{amidated_glycan_formula['H']}N{amidated_glycan_formula['N']}O{amidated_glycan_formula['O']}")
print(f"3. After Permethylation (x{num_methyl_sites}): C{final_glycan_formula['C']}H{final_glycan_formula['H']}N{final_glycan_formula['N']}O{final_glycan_formula['O']}")
print(f"\nThe mass of the final derivatized molecule (M) is {mass_M:.4f} Da.")
print(f"The observed ion is the singly-sodiated species [M+Na]+.")
print("\nFinal m/z Calculation:")
# Final output string as requested
final_eq = f"{mass_M:.4f} (Mass of M) + {mass_Na:.4f} (Mass of Na) = {mass_M_plus_Na:.4f} m/z"
print(final_eq)
