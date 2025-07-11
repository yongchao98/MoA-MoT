def calculate_molecular_formula():
    """
    Calculates the molecular formula of product A based on the reaction analysis.

    The reaction is an alkylation of a beta-keto ester, followed by saponification and decarboxylation.
    Product A is 2-benzyl-1-indanone.
    Its formula is calculated by taking the formula of the 1-indanone core (C9H8O),
    removing one hydrogen, and adding a benzyl group (C7H7).
    """

    # Define the atomic composition of the core structure and the substituent group
    indanone_core = {'C': 9, 'H': 8, 'O': 1}
    benzyl_group = {'C': 7, 'H': 7, 'O': 0}
    hydrogen_atom = {'C': 0, 'H': 1, 'O': 0}

    # Calculate the composition of the final product
    final_product = {}
    
    # Calculate carbons: Carbons from indanone core + Carbons from benzyl group
    final_product['C'] = indanone_core['C'] + benzyl_group['C']
    
    # Calculate hydrogens: Hydrogens from indanone core - 1 removed H + Hydrogens from benzyl group
    final_product['H'] = indanone_core['H'] - hydrogen_atom['H'] + benzyl_group['H']
    
    # Calculate oxygens: Oxygens from indanone core + Oxygens from benzyl group
    final_product['O'] = indanone_core['O'] + benzyl_group['O']

    # Print the step-by-step calculation
    print("Calculation for the molecular formula of product A (2-benzyl-1-indanone):")
    print(f"Carbon atoms = {indanone_core['C']} (from indanone core) + {benzyl_group['C']} (from benzyl group) = {final_product['C']}")
    print(f"Hydrogen atoms = {indanone_core['H']} (from indanone core) - {hydrogen_atom['H']} (removed) + {benzyl_group['H']} (from benzyl group) = {final_product['H']}")
    print(f"Oxygen atoms = {indanone_core['O']} (from indanone core) = {final_product['O']}")
    
    # Format and print the final molecular formula
    # Note: If an element count is 1, the number is usually omitted in the formula.
    c_str = f"C{final_product['C']}"
    h_str = f"H{final_product['H']}"
    o_str = "O" if final_product['O'] == 1 else f"O{final_product['O']}"

    molecular_formula = c_str + h_str + o_str
    
    print(f"\nThe molecular formula of compound A is: {molecular_formula}")


calculate_molecular_formula()