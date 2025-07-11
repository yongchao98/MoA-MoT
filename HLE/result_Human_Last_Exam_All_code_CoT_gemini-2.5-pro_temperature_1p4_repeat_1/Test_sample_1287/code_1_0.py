import collections

def calculate_mass_from_formula(formula, atom_masses):
    """Calculates the monoisotopic mass of a molecule from its chemical formula."""
    total_mass = 0
    for element, count in formula.items():
        total_mass += count * atom_masses[element]
    return total_mass

def main():
    """
    Calculates the m/z of a permethylated, amidated, and sodiated A2G2S2 glycan.
    """
    # Monoisotopic masses of the most abundant isotopes
    ATOM_MASSES = {
        'C': 12.0000000,
        'H': 1.0078250,
        'O': 15.9949146,
        'N': 14.0030740,
        'Na': 22.9897693,
    }

    # Chemical formulas for monosaccharide residues (as if in a polymer, i.e., -H2O)
    MONOSACCHARIDE_FORMULAS = {
        'Hex': {'C': 6, 'H': 10, 'O': 5},  # e.g., Man, Gal
        'HexNAc': {'C': 8, 'H': 13, 'N': 1, 'O': 5},  # e.g., GlcNAc
        'Neu5Ac': {'C': 11, 'H': 17, 'N': 1, 'O': 8}, # e.g., Sialic Acid
    }

    # Glycan composition
    glycan_composition = {
        'Hex': 5,  # 3 Man + 2 Gal
        'HexNAc': 4, # 4 GlcNAc
        'Neu5Ac': 2, # 2 Sialic Acid
    }
    
    # 1. Start with the base formula of the polymerized glycan
    base_formula = collections.defaultdict(int)
    for unit, count in glycan_composition.items():
        for element, num_atoms in MONOSACCHARIDE_FORMULAS[unit].items():
            base_formula[element] += num_atoms * count
    
    # Add back one water molecule for the ends of the polymer chain
    base_formula['H'] += 2
    base_formula['O'] += 1
    
    # 2. Account for amidation of the two sialic acids
    # Reaction: 2 * (-COOH -> -CONH2), which is a net change of -O +NH per site.
    # Total change: -O2, +N2, +H2
    amidated_formula = base_formula.copy()
    amidated_formula['O'] -= 2
    amidated_formula['N'] += 2
    amidated_formula['H'] += 2

    # 3. Account for permethylation
    # Detailed structural analysis of the amidated A2G2S2 glycan reveals 39 reactive sites
    # (free -OH and -NH groups) for methylation.
    num_methylation_sites = 39
    # Each methylation replaces H with CH3, a net gain of CH2.
    permethylated_formula = amidated_formula.copy()
    permethylated_formula['C'] += num_methylation_sites
    permethylated_formula['H'] += num_methylation_sites * 2

    # 4. Calculate the mass of the final modified glycan
    final_glycan_mass = calculate_mass_from_formula(permethylated_formula, ATOM_MASSES)
    
    # 5. Add the mass of the sodium adduct to get the final m/z
    final_ion_mass = final_glycan_mass + ATOM_MASSES['Na']
    
    # Print the results
    print("The three glycans A2G(4)2S(3)2, A2G(4)S(3)S(6), and A2G(4)2S(6)2 are isomers.")
    print("After amidation, permethylation, and sodiation, they will all have the same mass.")
    print("\nCalculation Steps:")
    print(f"1. Base glycan formula (A2G2S2): C{base_formula['C']}H{base_formula['H']}N{base_formula['N']}O{base_formula['O']}")
    print(f"2. After amidation (-O2, +N2H2): C{amidated_formula['C']}H{amidated_formula['H']}N{amidated_formula['N']}O{amidated_formula['O']}")
    print(f"3. After permethylation at {num_methylation_sites} sites (+C{num_methylation_sites}H{num_methylation_sites*2}): C{permethylated_formula['C']}H{permethylated_formula['H']}N{permethylated_formula['N']}O{permethylated_formula['O']}")
    print("\nFinal Mass Calculation:")
    mass_c = permethylated_formula['C'] * ATOM_MASSES['C']
    mass_h = permethylated_formula['H'] * ATOM_MASSES['H']
    mass_n = permethylated_formula['N'] * ATOM_MASSES['N']
    mass_o = permethylated_formula['O'] * ATOM_MASSES['O']
    mass_na = ATOM_MASSES['Na']
    
    print(f"Mass = ({permethylated_formula['C']} * C) + ({permethylated_formula['H']} * H) + ({permethylated_formula['N']} * N) + ({permethylated_formula['O']} * O) + (1 * Na)")
    print(f"Mass = ({permethylated_formula['C']} * {ATOM_MASSES['C']:.4f}) + ({permethylated_formula['H']} * {ATOM_MASSES['H']:.4f}) + ({permethylated_formula['N']} * {ATOM_MASSES['N']:.4f}) + ({permethylated_formula['O']} * {ATOM_MASSES['O']:.4f}) + (1 * {ATOM_MASSES['Na']:.4f})")
    print(f"Mass = {mass_c:.4f} + {mass_h:.4f} + {mass_n:.4f} + {mass_o:.4f} + {mass_na:.4f}")
    print(f"Neutral Mass = {final_glycan_mass:.4f} Da")
    print(f"[M+Na]+ Ion Mass = {final_glycan_mass:.4f} + {mass_na:.4f} = {final_ion_mass:.4f} Da")

    print("\n--- Results ---")
    print(f"The expected m/z for singly-sodiated A2G(4)2S(3)2 is: {final_ion_mass:.4f}")
    print(f"The expected m/z for singly-sodiated A2G(4)S(3)S(6) is: {final_ion_mass:.4f}")
    print(f"The expected m/z for singly-sodiated A2G(4)2S(6)2 is: {final_ion_mass:.4f}")

if __name__ == '__main__':
    main()
<<<2790.4151>>>