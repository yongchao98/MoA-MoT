import collections

def calculate_glycan_mass():
    """
    Calculates the m/z of a permethylated, amidated, and sodiated A2G2S2 glycan.
    """

    # Monoisotopic atomic masses
    ATOMIC_MASSES = {
        'C': 12.000000,
        'H': 1.007825,
        'O': 15.994915,
        'N': 14.003074,
        'Na': 22.989770,
    }

    # Chemical formulas of monosaccharide residues
    MONOSACCHARIDE_FORMULAS = {
        'Mannose': {'C': 6, 'H': 10, 'O': 5},        # C6H12O6 minus H2O
        'Galactose': {'C': 6, 'H': 10, 'O': 5},      # C6H12O6 minus H2O
        'GlcNAc': {'C': 8, 'H': 13, 'N': 1, 'O': 5}, # C8H15NO6 minus H2O
        'Neu5Ac': {'C': 11, 'H': 17, 'N': 1, 'O': 8},# C11H19NO9 minus H2O
    }
    
    # The A2G2S2 glycan is composed of a reducing-end GlcNAc and 10 other residues.
    # We add the formula for the final H2O back at the end.
    glycan_composition = {
        'Mannose': 3,
        'Galactose': 2,
        'GlcNAc': 4,
        'Neu5Ac': 2
    }

    print("Step 1: Calculating the chemical formula of the native A2G2S2 glycan.")
    # Summing up the atoms from all residues
    formula = collections.defaultdict(int)
    for name, count in glycan_composition.items():
        for atom, num in MONOSACCHARIDE_FORMULAS[name].items():
            formula[atom] += num * count
    
    # Add back the H2O for the non-reducing end glycan formula
    formula['H'] += 2
    formula['O'] += 1
    print(f"Native A2G2S2 Formula: C{formula['C']}H{formula['H']}N{formula['N']}O{formula['O']}")

    print("\nStep 2: Adjusting formula for DMT-MM amidation.")
    # Each of the 2 sialic acids has a -COOH converted to -CONH2.
    # Net change per reaction: -O, +N, +H. Total change: -2*O, +2*N, +2*H
    formula['O'] -= 2
    formula['N'] += 2
    formula['H'] += 2
    print(f"Amidated Glycan Formula: C{formula['C']}H{formula['H']}N{formula['N']}O{formula['O']}")
    
    print("\nStep 3: Calculating the number of sites for permethylation.")
    # Calculate active hydrogens on -OH and -NH groups that will be methylated.
    # -OH groups: Man(4), Gal(4), GlcNAc(3), Neu5Ac(4)
    oh_sites = (3 * 4) + (2 * 4) + (4 * 3) + (2 * 4)
    # -NH groups: GlcNAc(1), Neu5Ac acetyl(1), Neu5Ac primary amide(2)
    nh_sites = (4 * 1) + (2 * 1) + (2 * 2) 
    total_methylation_sites = oh_sites + nh_sites
    print(f"Total methylation sites (-OH and -NH): {total_methylation_sites}")

    print("\nStep 4: Adjusting formula for permethylation.")
    # Each methylation replaces H with CH3, a net addition of CH2.
    formula['C'] += total_methylation_sites
    formula['H'] += total_methylation_sites * 2
    final_formula = dict(sorted(formula.items()))
    print(f"Final derivatized glycan formula: C{final_formula['C']}H{final_formula['H']}N{final_formula['N']}O{final_formula['O']}")

    print("\nStep 5: Calculating the monoisotopic mass of the derivatized glycan (M).")
    mass_M = 0
    for atom, count in final_formula.items():
        mass_M += count * ATOMIC_MASSES[atom]
    print(f"Mass (M): {mass_M:.6f} Da")

    print("\nStep 6: Calculating the final m/z for the singly sodiated ion [M+Na]+.")
    mass_ion = mass_M + ATOMIC_MASSES['Na']
    charge = 1
    final_mz = mass_ion / charge
    
    print("\nSince all three glycans are isomers, they will have the same mass.")
    print("The theoretical m/z value for all three singly-sodiated species is:")
    
    # Final equation as requested
    c, h, n, o = final_formula['C'], final_formula['H'], final_formula['N'], final_formula['O']
    c_mass, h_mass, n_mass, o_mass, na_mass = ATOMIC_MASSES['C'], ATOMIC_MASSES['H'], ATOMIC_MASSES['N'], ATOMIC_MASSES['O'], ATOMIC_MASSES['Na']
    
    print(f"\nFinal Calculation:")
    print(f"({c} * {c_mass:.6f}) + ({h} * {h_mass:.6f}) + ({n} * {n_mass:.6f}) + ({o} * {o_mass:.6f}) + (1 * {na_mass:.6f}) = {final_mz:.6f}")
    
    return final_mz

# Run the calculation and store the final answer
result_mz = calculate_glycan_mass()
print(f"\n<<<A2G(4)2S(3)2: {result_mz:.4f}, A2G(4)S(3)S(6): {result_mz:.4f}, A2G(4)2S(6)2: {result_mz:.4f}>>>")
