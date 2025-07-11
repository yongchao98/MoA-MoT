def calculate_molecular_formula():
    """
    Calculates the molecular formula of the product through a three-step reaction.
    """
    # Step 1: Determine the molecular formula of the starting material.
    # Structure: 5-(Trifluoromethyl)-7-(p-methoxybenzyl)-7-azabicyclo[2.2.1]hept-5-en-2-one
    # Core (C6H5N): 6 C, 5 H, 1 N
    # C=O group: 1 O
    # CF3 group: 1 C, 3 F
    # PMB group (p-methoxybenzyl, -CH2-C6H4-OCH3): 8 C, 9 H, 1 O
    atoms = {
        'C': 6 + 1 + 8,  # Core + CF3 + PMB
        'H': 5 + 9,      # Core + PMB
        'F': 3,
        'N': 1,
        'O': 1 + 1       # Ketone + PMB ether
    }
    
    def format_formula(atom_dict):
        # Standard order: C, H, then alphabetical for others
        formula = ""
        if 'C' in atom_dict and atom_dict['C'] > 0:
            formula += f"C{atom_dict['C']}"
        if 'H' in atom_dict and atom_dict['H'] > 0:
            formula += f"H{atom_dict['H']}"
        
        sorted_keys = sorted([key for key in atom_dict if key not in ['C', 'H']])
        for key in sorted_keys:
            if atom_dict[key] > 0:
                formula += f"{key}{atom_dict[key]}"
        return formula

    print(f"Starting Material Formula: {format_formula(atoms)}")

    # Step 2: Reaction 1 - Deprotection of PMB group with CAN.
    # Remove PMB group (C8H9O) and add 1 H to the Nitrogen.
    atoms['C'] -= 8
    atoms['H'] -= 9
    atoms['O'] -= 1
    atoms['H'] += 1
    print(f"Intermediate 1 Formula (after deprotection): {format_formula(atoms)}")

    # Step 3: Reaction 2 - Hydrogenation with Pd/C, H2.
    # Reduce C=C double bond, adding 2 H atoms.
    atoms['H'] += 2
    print(f"Intermediate 2 Formula (after hydrogenation): {format_formula(atoms)}")

    # Step 4: Reaction 3 - Acid hydrolysis with 4 N HCl at 70 C.
    # This causes hydrolytic cleavage of a C-N bond in the strained bicyclic system.
    # This corresponds to the net addition of one water molecule (H2O).
    atoms['H'] += 2
    atoms['O'] += 1
    
    final_formula = format_formula(atoms)
    print(f"Final Product Formula (after hydrolysis): {final_formula}")
    
    # Print the final equation with each number explicitly.
    print("\nThe molecular formula of the product is determined as follows:")
    print(f"Number of Carbon atoms: {atoms['C']}")
    print(f"Number of Hydrogen atoms: {atoms['H']}")
    print(f"Number of Fluorine atoms: {atoms['F']}")
    print(f"Number of Nitrogen atoms: {atoms['N']}")
    print(f"Number of Oxygen atoms: {atoms['O']}")
    
calculate_molecular_formula()