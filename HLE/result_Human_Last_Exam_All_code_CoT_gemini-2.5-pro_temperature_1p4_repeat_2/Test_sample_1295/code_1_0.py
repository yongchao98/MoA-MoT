def solve_molecular_puzzle():
    """
    This function constructs and verifies a SMILES string based on a specific set of molecular properties.
    """

    # 1. Derived Molecular Formula and Proposed SMILES
    molecular_formula = "C8H18N6"
    smiles_representation = "NC(=N)C(C)(C)N=NC(C)(C)C(=N)N"

    # 2. Verification of key numerical properties, showing the "equation" as requested.
    print("--- Verification of Properties ---")

    # Property: Valence Electrons (Target: 80)
    num_C = 8
    num_H = 18
    num_N = 6
    valence_C = 4
    valence_H = 1
    valence_N = 5
    total_valence_electrons = (num_C * valence_C) + (num_H * valence_H) + (num_N * valence_N)
    print("\nValence Electron Calculation:")
    print(f"({num_C} C * {valence_C}) + ({num_H} H * {valence_H}) + ({num_N} N * {valence_N}) = {total_valence_electrons}")

    # Property: Molecular Weight (Target: 198.159)
    # Using monoisotopic masses for precision
    mass_C = 12.00000
    mass_H = 1.00783
    mass_N = 14.00307
    molecular_weight = (num_C * mass_C) + (num_H * mass_H) + (num_N * mass_N)
    print("\nMonoisotopic Molecular Weight Calculation:")
    print(f"({num_C} * {mass_C}) + ({num_H} * {mass_H:.5f}) + ({num_N} * {mass_N:.5f}) = {molecular_weight:.5f}")

    # 3. Final Answer Output
    print("\n--- Final Answer ---")
    print(f"The molecular formula is: {molecular_formula}")
    print("The constructed SMILES representation is:")
    print(smiles_representation)

# Execute the function to print the results
solve_molecular_puzzle()