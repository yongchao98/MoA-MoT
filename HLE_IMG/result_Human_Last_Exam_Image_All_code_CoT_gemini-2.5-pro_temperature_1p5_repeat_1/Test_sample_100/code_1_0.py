def solve_molecular_formula():
    """
    This function calculates the molecular formula of the product of a three-step reaction.
    """
    # Step 1: Define the composition of Intermediate 1.
    # The starting material is deprotected (PMB group C8H9O is removed, H is added)
    # The resulting structure (Intermediate 1) is a 6-(trifluoromethyl)-7-azabicyclo[2.2.1]hept-5-en-2-one.
    # Let's define its atomic composition.
    intermediate_1 = {'C': 7, 'H': 6, 'F': 3, 'N': 1, 'O': 1}
    print("The reaction sequence starts after the deprotection step.")
    print(f"Formula of Intermediate 1: C{intermediate_1['C']}H{intermediate_1['H']}F{intermediate_1['F']}N{intermediate_1['N']}O{intermediate_1['O']}")
    
    # Initialize the composition for tracking changes.
    current_composition = intermediate_1.copy()

    # Step 2: Hydrogenation (Pd/C, H2). This adds two hydrogen atoms.
    print("\nReaction 2: Hydrogenation adds 2 hydrogen atoms.")
    current_composition['H'] += 2
    print(f"Formula of Intermediate 2: C{current_composition['C']}H{current_composition['H']}F{current_composition['F']}N{current_composition['N']}O{current_composition['O']}")

    # Step 3: Acid Hydrolysis (HCl, H2O). This adds one water molecule.
    print("\nReaction 3: Lactam hydrolysis adds one water molecule (H2O).")
    current_composition['H'] += 2
    current_composition['O'] += 1

    # Final product composition
    final_composition = current_composition
    
    print("\nThe molecular formula of the final product is composed of:")
    print(f"Carbon atoms: {final_composition['C']}")
    print(f"Hydrogen atoms: {final_composition['H']}")
    print(f"Fluorine atoms: {final_composition['F']}")
    print(f"Nitrogen atoms: {final_composition['N']}")
    print(f"Oxygen atoms: {final_composition['O']}")
    
    # Construct the final formula string
    formula_str = (
        f"C{final_composition['C']}"
        f"H{final_composition['H']}"
        f"F{final_composition['F']}"
        f"N"  # N is 1, so no subscript
        f"O{final_composition['O']}"
    )
    print(f"\nThe molecular formula of the final product is: {formula_str}")

solve_molecular_formula()