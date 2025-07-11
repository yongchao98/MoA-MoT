def calculate_molecular_formula():
    """
    Calculates the molecular formula of the product of a three-step reaction.
    """

    # Step 1: Determine the formula of the starting material.
    # The structure is a 2-azabicyclo[2.2.2]oct-5-en-3-one derivative.
    # Core (deprotected, un-hydrogenated, C6-H): C7H9NO
    # Add CF3 substituent: C7H9NO - H + CF3 -> C8H8F3NO (This is Intermediate 1)
    # Add PMB group to N: C8H8F3NO - H + C8H9O (PMB group) -> C16H16F3NO2
    atoms = {'C': 16, 'H': 16, 'F': 3, 'N': 1, 'O': 2}

    print("Step-by-step calculation of the molecular formula:\n")

    def print_formula(atom_dict, name):
        # Helper function to print the formula in a standard format.
        order = ['C', 'H', 'F', 'N', 'O']
        formula_str = "".join([f"{el}{atom_dict.get(el, 0)}" for el in order if atom_dict.get(el, 0) > 0])
        print(f"Formula of {name}: {formula_str}")

    print("1. Starting Material:")
    print_formula(atoms, "Starting Material")
    print("-" * 20)

    # Step 2: Reaction 1 - Deprotection of PMB group using CAN.
    # This removes the p-methoxybenzyl group (C8H9O) and adds an H atom to the nitrogen.
    # Net change: C: -8, H: -8 (-9+1), O: -1
    atoms['C'] -= 8
    atoms['H'] -= 8
    atoms['O'] -= 1
    print("2. After PMB deprotection (forms Intermediate 1):")
    print_formula(atoms, "Intermediate 1")
    print("-" * 20)

    # Step 3: Reaction 2 - Hydrogenation of the C=C bond using Pd/C, H2.
    # This adds two hydrogen atoms.
    # Net change: H: +2
    atoms['H'] += 2
    print("3. After hydrogenation (forms Intermediate 2):")
    print_formula(atoms, "Intermediate 2")
    print("-" * 20)
    
    # Step 4: Reaction 3 - Hydrolysis of the lactam using aqueous HCl.
    # This adds one molecule of water (H2O) to break the amide bond.
    # Net change: H: +2, O: +1
    atoms['H'] += 2
    atoms['O'] += 1
    print("4. After lactam hydrolysis (forms Final Product):")
    print_formula(atoms, "Final Product")
    print("-" * 20)
    
    # Final Summary as requested.
    print("Final molecular formula composition:")
    print(f"Number of Carbon atoms (C): {atoms['C']}")
    print(f"Number of Hydrogen atoms (H): {atoms['H']}")
    print(f"Number of Fluorine atoms (F): {atoms['F']}")
    print(f"Number of Nitrogen atoms (N): {atoms['N']}")
    print(f"Number of Oxygen atoms (O): {atoms['O']}")

calculate_molecular_formula()