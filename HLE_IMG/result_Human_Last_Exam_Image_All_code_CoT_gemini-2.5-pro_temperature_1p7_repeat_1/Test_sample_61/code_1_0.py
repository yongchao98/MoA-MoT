def calculate_molecular_formula():
    """
    This script calculates the molecular formula for compound A.

    The reaction transforms methyl 1-oxoindane-2-carboxylate (1) and benzyl bromide (2)
    into 2-benzyl-1-indanone (A) through a sequence of alkylation, saponification,
    and decarboxylation.

    The final structure, 2-benzyl-1-indanone, can be thought of as a 1-indanone
    molecule where one hydrogen at the C2 position has been replaced by a benzyl group.

    We will calculate the molecular formula based on this final structure.
    - Molecular formula of 1-indanone: C9H8O
    - Molecular formula of a benzyl group: C7H7
    - An H atom is removed during substitution.
    """

    # Atom counts for the base molecule, 1-indanone
    indanone = {'C': 9, 'H': 8, 'O': 1}

    # Atom counts for the substituent, benzyl group
    benzyl_group = {'C': 7, 'H': 7, 'O': 0}

    # Atom counts for the replaced hydrogen atom
    replaced_h = {'C': 0, 'H': 1, 'O': 0}

    # Calculate the atom counts for the final product A
    final_C = indanone['C'] + benzyl_group['C']
    final_H = indanone['H'] + benzyl_group['H'] - replaced_h['H']
    final_O = indanone['O'] + benzyl_group['O']

    print("The final product A is 2-benzyl-1-indanone.")
    print("Its molecular formula is calculated by combining the atoms of 1-indanone and a benzyl group, and removing one hydrogen atom for the bond.")
    print("\nCalculation steps:")
    # Print the equation for each element
    print(f"Number of Carbon atoms (C): {indanone['C']} + {benzyl_group['C']} = {final_C}")
    print(f"Number of Hydrogen atoms (H): {indanone['H']} + {benzyl_group['H']} - {replaced_h['H']} = {final_H}")
    print(f"Number of Oxygen atoms (O): {indanone['O']} + {benzyl_group['O']} = {final_O}")

    # Format the final molecular formula string
    # We only include the 'O' if the count is > 0. Since it's 1, we don't show the number.
    formula = f"C{final_C}H{final_H}"
    if final_O > 0:
        if final_O > 1:
            formula += f"O{final_O}"
        else:
            formula += "O"
    
    print(f"\nThe molecular formula of compound A is {formula}.")


calculate_molecular_formula()
<<<C16H14O>>>