import collections

def identify_alkaloid():
    """
    Identifies the alkaloid compound from its structure, calculates its molecular
    weight, and prints the information.
    """
    # Name of the identified alkaloid
    compound_name = "Lupanine"

    # Chemical formula based on the structure
    chemical_formula = "C15H24N2O"

    # Atomic masses (g/mol)
    atomic_mass = {
        'C': 12.011,
        'H': 1.008,
        'N': 14.007,
        'O': 15.999
    }
    
    # Number of atoms for each element in the molecule
    atom_counts = {
        'C': 15,
        'H': 24,
        'N': 2,
        'O': 1
    }

    # Calculate molecular weight
    molecular_weight = (atom_counts['C'] * atomic_mass['C'] +
                        atom_counts['H'] * atomic_mass['H'] +
                        atom_counts['N'] * atomic_mass['N'] +
                        atom_counts['O'] * atomic_mass['O'])

    print(f"The name of the alkaloid is: {compound_name}")
    print(f"The chemical formula is: {chemical_formula}")
    print("-" * 30)
    print("Molecular Weight Calculation:")
    
    # Print the equation with each number as requested
    c_term = f"{atom_counts['C']} * {atomic_mass['C']}"
    h_term = f"{atom_counts['H']} * {atomic_mass['H']}"
    n_term = f"{atom_counts['N']} * {atomic_mass['N']}"
    o_term = f"{atom_counts['O']} * {atomic_mass['O']}"

    print(f"Equation: ({c_term}) + ({h_term}) + ({n_term}) + ({o_term})")
    
    print(f"Result: {molecular_weight:.3f} g/mol")

if __name__ == "__main__":
    identify_alkaloid()