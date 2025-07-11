def calculate_molecular_formula():
    """
    This function calculates the molecular formula of product A by summing the atoms
    of its constituent parts: the indanone core and the benzyl substituent.
    """
    # Initialize a dictionary to store the atom counts for product A
    atoms_A = {'C': 0, 'H': 0, 'O': 0}

    # Part 1: Define the atoms in the 1-indanone core (C9H8O) after removing
    # one hydrogen atom from C2 for substitution.
    indanone_core = {'C': 9, 'H': 8 - 1, 'O': 1}
    print("Step 1: Determine the atomic composition of the substituted 1-indanone core.")
    print(f"The 1-indanone core (C9H8O) loses one H for substitution, becoming: C={indanone_core['C']}, H={indanone_core['H']}, O={indanone_core['O']}\n")

    # Part 2: Define the atoms in the added benzyl group (C7H7).
    benzyl_group = {'C': 7, 'H': 7, 'O': 0}
    print("Step 2: Determine the atomic composition of the added substituent.")
    print(f"The benzyl group substituent is: C={benzyl_group['C']}, H={benzyl_group['H']}\n")

    # Step 3: Sum the atoms from the core and the substituent to get the final formula.
    atoms_A['C'] = indanone_core['C'] + benzyl_group['C']
    atoms_A['H'] = indanone_core['H'] + benzyl_group['H']
    atoms_A['O'] = indanone_core['O'] + benzyl_group['O']

    print("Step 3: Sum the atoms to find the molecular formula of the final product A.")
    print(f"Total Carbon (C): {indanone_core['C']} + {benzyl_group['C']} = {atoms_A['C']}")
    print(f"Total Hydrogen (H): {indanone_core['H']} + {benzyl_group['H']} = {atoms_A['H']}")
    print(f"Total Oxygen (O): {indanone_core['O']} + {benzyl_group['O']} = {atoms_A['O']}\n")

    # Print the final molecular formula in a standard format.
    final_formula = f"C{atoms_A['C']}H{atoms_A['H']}O{atoms_A['O']}"
    print(f"The final molecular formula of compound A is {final_formula}.")

# Run the calculation and print the result.
calculate_molecular_formula()