import collections

def calculate_glycan_mz():
    """
    Calculates the expected m/z for a specific series of glycan modifications.
    """
    # Step 1: Define constants and base composition
    # Monoisotopic atomic masses
    atomic_masses = {
        'C': 12.000000,
        'H': 1.007825,
        'O': 15.994915,
        'N': 14.003074,
        'Na': 22.989770
    }

    # Formulas for individual monosaccharides
    monosaccharide_formulas = {
        'Hex':      {'C': 6, 'H': 12, 'O': 6, 'N': 0},
        'HexNAc':   {'C': 8, 'H': 15, 'O': 6, 'N': 1},
        'Neu5Ac':   {'C': 11, 'H': 19, 'O': 9, 'N': 1}
    }
    
    # Glycan composition for A2G2S2
    glycan_composition = {
        'Hex': 5,
        'HexNAc': 4,
        'Neu5Ac': 2
    }
    
    print("--- Step 1: Calculating the formula of the native glycan (A2G2S2) ---")
    
    # Sum atoms from all monosaccharides
    base_formula = collections.defaultdict(int)
    num_residues = 0
    for residue, count in glycan_composition.items():
        num_residues += count
        for atom, atom_count in monosaccharide_formulas[residue].items():
            base_formula[atom] += atom_count * count

    # Subtract water for glycosidic bonds (n-1 bonds)
    num_bonds = num_residues - 1
    base_formula['H'] -= 2 * num_bonds
    base_formula['O'] -= 1 * num_bonds
    
    print(f"Glycan composition: {glycan_composition}")
    print(f"Total residues: {num_residues}, resulting in {num_bonds} glycosidic bonds.")
    print(f"Native glycan formula: C{base_formula['C']}H{base_formula['H']}N{base_formula['N']}O{base_formula['O']}\n")

    # Step 2: Apply amidation modification
    print("--- Step 2: Accounting for amidation of two sialic acids ---")
    # Reaction: 2x (-COOH -> -CONH2). Per reaction, change is -O, +N, +H.
    amidated_formula = base_formula.copy()
    amidated_formula['O'] -= 2
    amidated_formula['N'] += 2
    amidated_formula['H'] += 2
    print("Reaction converts 2 carboxyls to primary amides.")
    print(f"Formula change: -2 O, +2 N, +2 H")
    print(f"Amidated glycan formula: C{amidated_formula['C']}H{amidated_formula['H']}N{amidated_formula['N']}O{amidated_formula['O']}\n")

    # Step 3: Apply permethylation modification
    print("--- Step 3: Accounting for permethylation ---")
    # Count methylation sites in the amidated glycan structure.
    # This includes all free -OH groups and all N-H protons.
    # Detailed count: 29 (-OH) + 10 (N-H) = 39 sites.
    # (4 GlcNAc N-H + 2 Neu5Ac N-H + 4 N-H from the two new amides = 10 N-H)
    num_methylations = 39
    
    # Permethylation replaces H with CH3, a net addition of CH2.
    final_formula = amidated_formula.copy()
    final_formula['C'] += num_methylations
    final_formula['H'] += 2 * num_methylations
    
    print(f"Number of methylation sites (free -OH and N-H groups): {num_methylations}")
    print(f"Formula change: +{num_methylations} C, +{2 * num_methylations} H")
    print(f"Final modified glycan formula: C{final_formula['C']}H{final_formula['H']}N{final_formula['N']}O{final_formula['O']}\n")
    
    # Step 4: Calculate final m/z
    print("--- Step 4: Calculating the final m/z for the [M+Na]+ ion ---")
    
    # Helper function to calculate mass from a formula dictionary
    def calculate_mass(formula, masses):
        mass = sum(formula[atom] * masses[atom] for atom in formula)
        return mass

    mass_M = calculate_mass(final_formula, atomic_masses)
    mass_Na = atomic_masses['Na']
    final_mz = mass_M + mass_Na

    print(f"The isomers A2G(4)2S(3)2, A2G(4)S(3)S(6), and A2G(4)2S(6)2 have the same mass.")
    print("Final m/z is calculated as: Mass(Final Glycan) + Mass(Na)")
    print(f"Final m/z = {mass_M:.4f} + {mass_Na:.4f}")
    print(f"\nThe expected m/z for all three singly-sodiated glycans is: {final_mz:.4f}")
    
    return final_mz

expected_mz = calculate_glycan_mz()
# The final answer is wrapped according to the required format.
# The value is formatted to 4 decimal places, which is standard for high-resolution MS.
print(f"\n<<<Since A2G(4)2S(3)2, A2G(4)S(3)S(6), and A2G(4)2S(6)2 are isomers, they all have an expected m/z of {expected_mz:.4f} for the singly-sodiated ion.>>>")