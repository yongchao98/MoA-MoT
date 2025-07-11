def calculate_glycan_mass():
    """
    Calculates the m/z of specific permethylated, amidated, and sodiated N-glycans.
    """

    # --- Step 1: Define constants and composition ---
    
    # Monoisotopic atomic masses
    atomic_mass = {
        'C': 12.000000,
        'H': 1.007825,
        'O': 15.994915,
        'N': 14.003074,
        'Na': 22.989770
    }

    # Base chemical formulas for each monomer
    monomer_formula = {
        'Hex':    {'C': 6, 'H': 12, 'O': 6, 'N': 0}, # e.g., Galactose, Mannose
        'HexNAc': {'C': 8, 'H': 15, 'O': 6, 'N': 1}, # e.g., GlcNAc
        'NeuAc':  {'C': 11, 'H': 19, 'O': 9, 'N': 1} # Sialic Acid
    }

    # Composition of the A2G2S2 glycan
    glycan_composition = {
        'Hex': 5,     # 2 Galactose + 3 Mannose
        'HexNAc': 4,  # 4 GlcNAc
        'NeuAc': 2     # 2 Sialic Acids
    }
    
    glycan_names = ["A2G(4)2S(3)2", "A2G(4)S(3)S(6)", "A2G(4)2S(6)2"]

    # --- Step 2: Calculate base molecular formula of the un-modified glycan ---

    # Sum formulas of all monomers
    base_formula = {'C': 0, 'H': 0, 'O': 0, 'N': 0}
    num_residues = 0
    for key, count in glycan_composition.items():
        num_residues += count
        for atom, atom_count in monomer_formula[key].items():
            base_formula[atom] += atom_count * count

    # Subtract water for each glycosidic bond (N-1 bonds for N residues)
    num_bonds = num_residues - 1
    base_formula['H'] -= 2 * num_bonds
    base_formula['O'] -= 1 * num_bonds

    # --- Step 3: Apply chemical modification formulas ---
    
    # a) Amidation of the 2 sialic acid carboxyl groups (COOH -> CONH2)
    # Net change per reaction: -O, +N, +H
    amidation_modifications = glycan_composition['NeuAc']
    modified_formula = base_formula.copy()
    modified_formula['O'] -= 1 * amidation_modifications
    modified_formula['N'] += 1 * amidation_modifications
    modified_formula['H'] += 1 * amidation_modifications

    # b) Permethylation of all free -OH and -NH groups
    # The number of sites must be calculated.
    # Labile H in un-modified glycan:
    # - OH groups: (5 Hex * 5 OH) + (4 HexNAc * 3 OH) + (2 NeuAc * 5 OH) - 2*num_bonds = 25+12+10 - 2*10 = 27
    # - NH groups (from N-acetyl): (4 HexNAc * 1 NH) + (2 NeuAc * 1 NH) = 6
    # - COOH groups: (2 NeuAc * 1 COOH) = 2
    # Total labile H before amidation = 27 + 6 + 2 = 35
    # After amidation, 2 COOH (2 H) are replaced by 2 CONH2 (4 H). Net gain of 2 labile H.
    # Total methylation sites = 35 - 2 + 4 = 37
    num_methylation_sites = 37
    
    # Add a CH2 group for each methylation (-H replaced by -CH3)
    modified_formula['C'] += 1 * num_methylation_sites
    modified_formula['H'] += 2 * num_methylation_sites

    final_formula = modified_formula
    
    # --- Step 4 & 5: Calculate the final mass of the sodiated ion [M+Na]+ ---
    
    neutral_mass = 0
    for atom, count in final_formula.items():
        neutral_mass += count * atomic_mass[atom]
        
    sodiated_ion_mass = neutral_mass + atomic_mass['Na']

    # --- Step 6: Print the results ---

    print(f"The final molecular formula of the modified glycan is C{final_formula['C']}H{final_formula['H']}N{final_formula['N']}O{final_formula['O']}.\n")
    print("The theoretical m/z of the singly-sodiated ion [M+Na]+ is calculated as follows:\n")
    
    mass_eq_str = (
        f"Mass = {final_formula['C']}*C + {final_formula['H']}*H + {final_formula['N']}*N + "
        f"{final_formula['O']}*O + Na\n"
        f"     = {final_formula['C']}*{atomic_mass['C']} + {final_formula['H']}*{atomic_mass['H']:.6f} + "
        f"{final_formula['N']}*{atomic_mass['N']:.6f} + {final_formula['O']}*{atomic_mass['O']:.6f} + {atomic_mass['Na']:.6f}"
    )
    print(mass_eq_str)
    print(f"     = {sodiated_ion_mass:.4f}\n")
    
    print("Since sialic acid linkage isomers have the same elemental composition, the expected masses are:")
    for name in glycan_names:
        print(f"The expected m/z for {name} is {sodiated_ion_mass:.4f}")

    return sodiated_ion_mass

if __name__ == '__main__':
    final_mass = calculate_glycan_mass()
    # The final answer is wrapped for evaluation purposes
    print(f"\n<<<{final_mass:.4f}>>>")
