def solve_molecule_design():
    """
    This script presents the solution for the molecule design problem and
    verifies its molecular weight based on the derived formula.
    """
    # The SMILES representation of the designed molecule
    smiles_string = "O1CCN(CC1)CCOCCN2CCOCC2"

    # Molecular formula derived from the constraints
    formula = "C12H24N2O3"
    atom_counts = {'C': 12, 'H': 24, 'N': 2, 'O': 3}

    # Exact atomic masses for calculation
    atomic_mass = {
        'C': 12.00000,
        'H': 1.007825,
        'N': 14.003074,
        'O': 15.994915
    }

    print(f"Proposed SMILES string: {smiles_string}")
    print(f"Corresponding Molecular Formula: {formula}\n")
    print("Verifying the molecular weight against the target of 244.179:")

    # Calculate contribution of each element to the total molecular weight
    c_mass = atom_counts['C'] * atomic_mass['C']
    h_mass = atom_counts['H'] * atomic_mass['H']
    n_mass = atom_counts['N'] * atomic_mass['N']
    o_mass = atom_counts['O'] * atomic_mass['O']

    # Calculate the total molecular weight
    total_mw = c_mass + h_mass + n_mass + o_mass

    # Print the equation as requested
    print("\nMolecular Weight Calculation Equation:")
    print(f"({atom_counts['C']} * {atomic_mass['C']}) + ({atom_counts['H']} * {atomic_mass['H']}) + ({atom_counts['N']} * {atomic_mass['N']}) + ({atom_counts['O']} * {atomic_mass['O']})")
    print(f"= {c_mass:.5f} + {h_mass:.5f} + {n_mass:.5f} + {o_mass:.5f}")
    print(f"= {total_mw:.5f}\n")

    print(f"The calculated molecular weight is {total_mw:.5f}, which matches the required value.")
    print("\nThe final SMILES representation that fulfills all conditions is:")
    print(smiles_string)

solve_molecule_design()
<<<O1CCN(CC1)CCOCCN2CCOCC2>>>