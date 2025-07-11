def solve_chemistry_problem():
    """
    This function identifies the final product of the described chemical synthesis.
    The final product, Compound 4, is Phthalic Acid.
    This script calculates and prints its properties.
    """

    # --- Properties of Phthalic Acid ---
    compound_name = "Phthalic acid"
    
    # Molecular Formula: C8H6O4
    atom_counts = {
        'Carbon': 8,
        'Hydrogen': 6,
        'Oxygen': 4
    }
    
    molecular_formula = f"C{atom_counts['Carbon']}H{atom_counts['Hydrogen']}O{atom_counts['Oxygen']}"

    # Atomic weights (g/mol)
    atomic_weights = {
        'C': 12.011,
        'H': 1.008,
        'O': 15.999
    }

    # Calculate Molecular Weight
    molecular_weight = (atom_counts['Carbon'] * atomic_weights['C'] +
                        atom_counts['Hydrogen'] * atomic_weights['H'] +
                        atom_counts['Oxygen'] * atomic_weights['O'])

    # --- Print the final answer ---
    print(f"The final product, Compound 4, is: {compound_name}")
    print(f"Molecular Formula: {molecular_formula}")
    print("\n--- Atoms in the final molecule (equation) ---")
    print(f"Number of Carbon atoms: {atom_counts['Carbon']}")
    print(f"Number of Hydrogen atoms: {atom_counts['Hydrogen']}")
    print(f"Number of Oxygen atoms: {atom_counts['Oxygen']}")
    print(f"\nCalculated Molecular Weight: {molecular_weight:.3f} g/mol")

solve_chemistry_problem()