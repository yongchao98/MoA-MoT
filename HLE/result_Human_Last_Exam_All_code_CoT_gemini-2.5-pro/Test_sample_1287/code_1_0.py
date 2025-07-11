def calculate_glycan_mass():
    """
    Calculates the m/z of a derivatized biantennary disialylated N-glycan.

    The glycan undergoes two reactions:
    1. Amidation of both sialic acids.
    2. Permethylation of all free hydroxyl and amide groups.

    The final ion observed is the singly-sodiated species [M+Na]+.
    """
    # Monoisotopic atomic masses
    atomic_masses = {
        'C': 12.000000,
        'H': 1.007825,
        'O': 15.994915,
        'N': 14.003074,
        'Na': 22.989770
    }

    # Step 1: Define the initial glycan composition and calculate its mass.
    # The base glycan is a biantennary glycan with composition:
    # 3x Mannose, 4x N-acetylglucosamine, 2x Galactose, 2x N-acetylneuraminic acid (sialic acid)
    # Total composition: Man(3)GlcNAc(4)Gal(2)Neu5Ac(2)
    # The molecular formula for the free glycan (with 10 glycosidic bonds) is C84 H138 N6 O62.
    initial_formula = {'C': 84, 'H': 138, 'N': 6, 'O': 62}
    
    initial_mass = sum(atomic_masses[atom] * count for atom, count in initial_formula.items())

    print("--- Step 1: Initial Glycan Mass ---")
    print(f"The three glycans are isomers with the same composition: Man(3)GlcNAc(4)Gal(2)Neu5Ac(2)")
    print(f"Initial molecular formula: C84H138N6O62")
    print(f"Initial neutral mass: {initial_mass:.4f} Da\n")

    # Step 2: Calculate the mass change from amidation.
    # The reaction converts 2 carboxylic acids (-COOH) to primary amides (-CONH2).
    # Per reaction, the change is -O, +N, +H. For two reactions: -2O, +2N, +2H.
    amidation_formula_change = {'H': 2, 'N': 2, 'O': -2}
    amidated_formula = {atom: initial_formula.get(atom, 0) + amidation_formula_change.get(atom, 0) for atom in set(initial_formula) | set(amidation_formula_change)}
    amidated_mass = sum(atomic_masses[atom] * count for atom, count in amidated_formula.items())
    
    print("--- Step 2: Amidation ---")
    print("Reaction: 2x -COOH -> 2x -CONH2")
    print(f"Formula after amidation: C{amidated_formula['C']}H{amidated_formula['H']}N{amidated_formula['N']}O{amidated_formula['O']}")
    print(f"Neutral mass after amidation: {amidated_mass:.4f} Da\n")
    
    # Step 3: Calculate the mass change from permethylation.
    # A standard disialylated biantennary glycan has 37 sites for methylation (-OH and -NH groups),
    # assuming ring opening at the reducing end.
    # Each methylation replaces an H with a CH3, a net addition of CH2.
    num_methylation_sites = 37
    permethylation_formula_change = {'C': num_methylation_sites, 'H': 2 * num_methylation_sites}
    final_neutral_formula = {atom: amidated_formula.get(atom, 0) + permethylation_formula_change.get(atom, 0) for atom in set(amidated_formula) | set(permethylation_formula_change)}
    final_neutral_mass = sum(atomic_masses[atom] * count for atom, count in final_neutral_formula.items())
    
    print("--- Step 3: Permethylation ---")
    print(f"Reaction: Methylation at {num_methylation_sites} sites (all -OH and -NH groups)")
    print(f"Formula after permethylation: C{final_neutral_formula['C']}H{final_neutral_formula['H']}N{final_neutral_formula['N']}O{final_neutral_formula['O']}")
    print(f"Final neutral mass (M): {final_neutral_mass:.4f} Da\n")

    # Step 4: Calculate the final m/z for the sodiated ion [M+Na]+.
    mass_na = atomic_masses['Na']
    final_mz = final_neutral_mass + mass_na

    print("--- Step 4: Final m/z Calculation ---")
    print("The mass is the same for A2G(4)2S(3)2, A2G(4)S(3)S(6), and A2G(4)2S(6)2 as they are isomers.")
    print("The final observed ion is the singly-sodiated species [M+Na]+.")
    print("\nFinal Equation:")
    print(f"Mass (M) + Mass (Na) = m/z ([M+Na]+)")
    print(f"{final_neutral_mass:.4f} + {mass_na:.4f} = {final_mz:.4f}")
    
    return final_mz

# Run the calculation and store the final answer
expected_mz = calculate_glycan_mass()

# The final answer is wrapped as requested
print(f"\n<<<The expected m/z for all three glycans is {expected_mz:.4f}>>>")