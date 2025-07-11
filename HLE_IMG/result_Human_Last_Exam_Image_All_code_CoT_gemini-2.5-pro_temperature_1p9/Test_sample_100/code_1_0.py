def calculate_molecular_formula():
    """
    Calculates the molecular formula of the final product step-by-step.
    """
    # Step 1: Define the atomic composition of the starting material
    # Formula: C15H14F3NO2
    atoms = {'C': 15, 'H': 14, 'F': 3, 'N': 1, 'O': 2}
    print("Step 0: Starting Material Formula")
    print_formula(atoms)

    # Step 2: Reaction 1 - PMB deprotection with CAN
    # The PMB group (p-methoxybenzyl, C8H9O) is removed, and an H atom is added.
    # Change: -C8, -H9, -O1, +H1  (Net: -C8, -H8, -O1)
    atoms['C'] -= 8
    atoms['H'] -= 9
    atoms['O'] -= 1
    atoms['H'] += 1
    print("\nStep 1: After PMB deprotection (Intermediate 1)")
    print_formula(atoms)

    # Step 3: Reaction 2 - Hydrogenation with Pd/C, H2
    # The C=C double bond is reduced, adding two H atoms.
    # Change: +H2
    atoms['H'] += 2
    print("\nStep 2: After hydrogenation (Intermediate 2)")
    print_formula(atoms)

    # Step 4: Reaction 3 - Amide hydrolysis with HCl
    # The lactam is hydrolyzed, consuming one molecule of water.
    # Change: +H2O
    atoms['H'] += 2
    atoms['O'] += 1
    print("\nStep 3: After amide hydrolysis (Final Product)")
    print_formula(atoms)

def print_formula(atom_dict):
    """Helper function to print the molecular formula."""
    formula = f"C{atom_dict['C']}H{atom_dict['H']}F{atom_dict['F']}NO{atom_dict['O']}"
    print(f"The molecular formula is: {formula}")

calculate_molecular_formula()