def solve_molecule_puzzle():
    """
    This function provides the solution to the molecule design puzzle.
    It prints the final SMILES string and the detailed molecular weight calculation.
    """

    # The final SMILES string for the molecule that satisfies all constraints.
    final_smiles = "CCN=C(C)Cn1cnc(c1)c2ccc(O)cc2"

    # Molecular formula: C14H17N3O
    c_count = 14
    h_count = 17
    n_count = 3
    o_count = 1

    # Exact atomic masses
    c_mass = 12.00000
    h_mass = 1.007825
    n_mass = 14.003074
    o_mass = 15.994915

    # Calculate the total molecular weight
    total_mw = (c_count * c_mass) + (h_count * h_mass) + (n_count * n_mass) + (o_count * o_mass)

    print(f"Final SMILES: {final_smiles}")
    print("\nMolecular Weight Calculation:")
    
    # Print the equation with each number, as requested.
    equation = (
        f"{c_count} * {c_mass:.5f} (C) + "
        f"{h_count} * {h_mass:.6f} (H) + "
        f"{n_count} * {n_mass:.6f} (N) + "
        f"{o_count} * {o_mass:.6f} (O) = "
        f"{total_mw:.5f}"
    )
    print(equation)

if __name__ == "__main__":
    solve_molecule_puzzle()
<<<CCN=C(C)Cn1cnc(c1)c2ccc(O)cc2>>>