def final_reaction_equation():
    """
    This function stores and prints the stoichiometric numbers from the final step of the proposed synthesis.
    The reaction is the conversion of the precursor to the target molecule:
    C₃₁H₂₂O₂ -> C₃₀H₁₄ + CO + 4H₂
    """
    precursor_carbons = 31
    precursor_hydrogens = 22
    precursor_oxygens = 2
    
    product_carbons = 30
    product_hydrogens = 14
    
    co_carbons = 1
    co_oxygens = 1
    
    h2_hydrogens_per_molecule = 2
    h2_molecules = 4
    
    print("The numbers in the final chemical equation (C₃₁H₂₂O₂ -> C₃₀H₁₄ + CO + 4H₂) are:")
    # Printing each number from the equation's molecules
    print(precursor_carbons, precursor_hydrogens, precursor_oxygens, "->", 
          product_carbons, product_hydrogens, "+", 
          co_carbons, co_oxygens, "+", 
          h2_molecules, "*", h2_hydrogens_per_molecule)

final_reaction_equation()