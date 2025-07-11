import collections

def calculate_mass_from_formula(formula, monoisotopic_masses):
    """Calculates the monoisotopic mass from a chemical formula dictionary."""
    mass = 0.0
    for element, count in formula.items():
        mass += count * monoisotopic_masses[element]
    return mass

def main():
    """
    Calculates the m/z of a permethylated, amidated, and sodiated A2G2S2 glycan.
    """
    # 1. Define constants
    MONOISOTOPIC_MASSES = {
        'C': 12.000000,
        'H': 1.007825,
        'O': 15.994915,
        'N': 14.003074,
        'Na': 22.989770,
    }

    # 2. Define monomer chemical formulas
    MONOMER_FORMULAS = {
        'Hex': {'C': 6, 'H': 12, 'O': 6},        # Mannose, Galactose
        'HexNAc': {'C': 8, 'H': 15, 'N': 1, 'O': 6}, # N-acetylglucosamine
        'Neu5Ac': {'C': 11, 'H': 19, 'N': 1, 'O': 9} # N-acetylneuraminic acid
    }
    
    # Define composition of A2G2S2 glycan
    GLYCAN_COMPOSITION = {
        'Hex': 5,     # 3 Man + 2 Gal
        'HexNAc': 4,
        'Neu5Ac': 2
    }
    
    num_residues = sum(GLYCAN_COMPOSITION.values())
    num_glycosidic_bonds = num_residues - 1

    # 3. Calculate chemical formula of the native, free glycan
    native_formula = collections.defaultdict(int)
    for monomer, count in GLYCAN_COMPOSITION.items():
        for element, num_atoms in MONOMER_FORMULAS[monomer].items():
            native_formula[element] += count * num_atoms
    
    # Account for water loss from forming glycosidic bonds
    native_formula['H'] -= 2 * num_glycosidic_bonds
    native_formula['O'] -= 1 * num_glycosidic_bonds

    mass_native = calculate_mass_from_formula(native_formula, MONOISOTOPIC_MASSES)
    
    # 4. Apply amidation modification to the formula
    # Reaction: -COOH -> -CONH2 per sialic acid. Net change: -O, +N, +H
    # Two sialic acids are amidated.
    amidated_formula = native_formula.copy()
    amidated_formula['O'] -= 2 * 1
    amidated_formula['N'] += 2 * 1
    amidated_formula['H'] += 2 * 1
    
    mass_amidated = calculate_mass_from_formula(amidated_formula, MONOISOTOPIC_MASSES)
    
    # 5. Apply permethylation modification to the formula
    # Count methylation sites on the amidated glycan:
    # Sites on free monomers: Hex(5 OH), HexNAc(4 OH), Amidated-Neu5Ac(5 OH + 2 NH)
    sites_hex = 5
    sites_hexnac = 4
    sites_amidated_neu5ac = 7 # 5 -OH groups + 2 N-H from the new amide
    
    total_sites_on_monomers = (GLYCAN_COMPOSITION['Hex'] * sites_hex +
                               GLYCAN_COMPOSITION['HexNAc'] * sites_hexnac +
                               GLYCAN_COMPOSITION['Neu5Ac'] * sites_amidated_neu5ac)
                               
    num_methylations = total_sites_on_monomers - (2 * num_glycosidic_bonds)
    
    # Modify formula for permethylation: replace H with CH3 (net change +CH2)
    # The H is replaced, so we add 3*H and subtract 1*H = add 2*H
    final_formula = amidated_formula.copy()
    final_formula['C'] += num_methylations * 1
    final_formula['H'] += num_methylations * 2

    mass_final_molecule = calculate_mass_from_formula(final_formula, MONOISOTOPIC_MASSES)

    # 6. Calculate final m/z for the [M+Na]+ ion
    mz_value = mass_final_molecule + MONOISOTOPIC_MASSES['Na']

    # --- Output Results ---
    print("The three glycans A2G(4)2S(3)2, A2G(4)S(3)S(6), and A2G(4)2S(6)2 are isomers.")
    print("Isomers have the same chemical formula and therefore the same mass.")
    print("The expected m/z value will be identical for all three.\n")
    print("--- Calculation Steps ---")
    
    formula_str_native = "".join([f"{el}{count}" for el, count in sorted(native_formula.items())])
    print(f"1. Native A2G2S2 Glycan Formula: {formula_str_native}")
    print(f"   Mass of Native Glycan: {mass_native:.4f} Da\n")
    
    print(f"2. After Amidation of 2 Sialic Acids (-COOH -> -CONH2):")
    formula_str_amidated = "".join([f"{el}{count}" for el, count in sorted(amidated_formula.items())])
    print(f"   Formula becomes: {formula_str_amidated}")
    print(f"   Mass of Amidated Glycan: {mass_amidated:.4f} Da\n")

    print(f"3. After Permethylation:")
    print(f"   Number of methylation sites calculated: {num_methylations}")
    formula_str_final = "".join([f"{el}{count}" for el, count in sorted(final_formula.items())])
    print(f"   Final permethylated, amidated formula: {formula_str_final}")
    print(f"   Mass of Final Molecule (M): {mass_final_molecule:.4f} Da\n")
    
    print("4. Final m/z Calculation for Singly-Sodiated Ion [M+Na]+:")
    print("   m/z = Mass(M) + Mass(Na+)")
    print(f"   m/z = {mass_final_molecule:.4f} + {MONOISOTOPIC_MASSES['Na']:.4f}")
    print(f"   Final Observed m/z = {mz_value:.4f}\n")
    
    print("Therefore, the mass you should observe is the same for all three isomers.")

if __name__ == '__main__':
    main()