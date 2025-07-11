def calculate_molecular_formula():
    """
    Calculates the molecular formula of product A based on the reaction provided.
    The reaction is a phase-transfer catalyzed alkylation of a beta-keto ester,
    followed by saponification and decarboxylation.
    """

    # Step 1: Define the molecular formula of reactant 1 (methyl 1-oxo-2,3-dihydro-1H-indene-2-carboxylate)
    # It consists of an indanone core (C9H8O) where one H at C2 is replaced by a methoxycarbonyl group (COOCH3).
    # Indanone: C9H8O. Replacing H with COOMe means C9H7O + C2H3O2 = C11H10O3
    c1_formula = {'C': 11, 'H': 10, 'O': 3}
    print(f"Molecular formula of compound 1: C{c1_formula['C']}H{c1_formula['H']}O{c1_formula['O']}")

    # The reaction involves two main transformations to get to product A.

    # Transformation 1: Alkylation
    # The acidic hydrogen at the alpha-position of compound 1 is replaced by a benzyl group (C7H7) from compound 2.
    # We add the atoms of a benzyl group and subtract one hydrogen atom.
    print("\nStep 1: Alkylation with a benzyl group (C7H7)")
    intermediate_formula = c1_formula.copy()
    benzyl_group = {'C': 7, 'H': 7}
    intermediate_formula['C'] += benzyl_group['C']
    intermediate_formula['H'] += benzyl_group['H'] - 1
    print(f"  - Add benzyl group: +{benzyl_group['C']} C, +{benzyl_group['H']} H")
    print("  - Remove 1 alpha-H: -1 H")
    print(f"Formula of alkylated intermediate: C{intermediate_formula['C']}H{intermediate_formula['H']}O{intermediate_formula['O']}")

    # Transformation 2: Saponification and Decarboxylation
    # The methoxycarbonyl group (-COOMe, or -C2H3O2) is removed and replaced with a hydrogen atom.
    print("\nStep 2: Saponification and Decarboxylation")
    final_formula = intermediate_formula.copy()
    methoxycarbonyl_group = {'C': 2, 'H': 3, 'O': 2}
    final_formula['C'] -= methoxycarbonyl_group['C']
    final_formula['H'] -= methoxycarbonyl_group['H']
    final_formula['O'] -= methoxycarbonyl_group['O']
    final_formula['H'] += 1 # Add one hydrogen back
    print(f"  - Remove -COOMe group: -{methoxycarbonyl_group['C']} C, -{methoxycarbonyl_group['H']} H, -{methoxycarbonyl_group['O']} O")
    print("  - Add 1 hydrogen back: +1 H")

    # Print the final molecular formula of compound A
    print("\n--------------------------------------------------")
    print("Final Molecular Formula of Compound A:")
    print(f"Carbon atoms = {final_formula['C']}")
    print(f"Hydrogen atoms = {final_formula['H']}")
    print(f"Oxygen atoms = {final_formula['O']}")
    print(f"The molecular formula of compound A is C{final_formula['C']}H{final_formula['H']}O{final_formula['O']}.")
    print("--------------------------------------------------")

calculate_molecular_formula()