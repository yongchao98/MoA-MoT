def identify_compound_C():
    """
    This function identifies and describes the chemical compound C.
    """
    # Name of the final compound C
    name_C = "bis(2,4,6-trihydroxyphenyl)(4-diethylamino-2,6-dihydroxyphenyl)methylium"

    # The structure of compound C in SMILES format
    smiles_C = "[C+](c1c(O)cc(O)cc1O)(c1c(O)cc(O)cc1O)(c1c(O)cc(N(CC)CC)cc1O)"

    # Atom counts for the molecular formula calculation ("the numbers in the final equation")
    carbon_atoms = 23
    hydrogen_atoms = 24
    nitrogen_atoms = 1
    oxygen_atoms = 8
    charge = "+1"

    # Constructing the molecular formula string
    molecular_formula = f"C{carbon_atoms}H{hydrogen_atoms}N{nitrogen_atoms}O{oxygen_atoms}"

    # Printing the results
    print("Identity of Final Compound C:")
    print("-" * 30)
    print(f"Name: {name_C}")
    print(f"SMILES String: {smiles_C}")
    print("\n" + "-"*30)
    print("Molecular Formula Calculation:")
    print(f"Number of Carbon atoms = {carbon_atoms}")
    print(f"Number of Hydrogen atoms = {hydrogen_atoms}")
    print(f"Number of Nitrogen atoms = {nitrogen_atoms}")
    print(f"Number of Oxygen atoms = {oxygen_atoms}")
    print(f"Overall Charge = {charge}")
    print("-" * 30)
    print(f"Final Molecular Formula: [{molecular_formula}]{charge}")
    print("-" * 30)


# Run the function to display the answer
identify_compound_C()
