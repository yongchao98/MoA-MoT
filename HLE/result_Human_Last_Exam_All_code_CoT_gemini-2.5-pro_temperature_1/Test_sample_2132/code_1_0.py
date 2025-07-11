def calculate_bottcher_complexity():
    """
    Calculates the Böttcher Molecular Complexity for the product of the
    Favorskii rearrangement of 2-chlorocyclohexanone.
    """
    # Step 1: Identify the product of the reaction.
    # The Favorskii rearrangement of 2-chlorocyclohexanone results in the
    # ring-contracted product, cyclopentanecarboxylic acid (C6H10O2).
    product_name = "cyclopentanecarboxylic acid"
    product_formula = "C6H10O2"

    print(f"The reaction is the Favorskii rearrangement of 2-chlorocyclohexanone.")
    print(f"The product of this reaction is {product_name}, with the formula {product_formula}.\n")

    # Step 2: Define the components for the Böttcher Molecular Complexity formula.
    # Formula: MC = (N_atoms * N_bonds) / N_parts

    # N_atoms: Total number of atoms in the molecule.
    # C=6, H=10, O=2
    n_atoms = 6 + 10 + 2

    # N_bonds: Sum of bond orders in the molecule.
    # Cyclopentanecarboxylic acid has:
    # - 6 C-C single bonds
    # - 9 C-H single bonds
    # - 1 C=O double bond (bond order = 2)
    # - 1 C-O single bond
    # - 1 O-H single bond
    # Total = 6 + 9 + 2 + 1 + 1 = 19
    n_bonds = 19

    # N_parts: Number of disconnected parts of the molecule. For a single molecule, this is 1.
    n_parts = 1

    print("Calculating the Böttcher Molecular Complexity:")
    print(f"Formula: MC = (Number of Atoms * Number of Bonds) / Number of Parts")
    print(f"Number of Atoms (N_atoms) = {n_atoms}")
    print(f"Number of Bonds (N_bonds) = {n_bonds}")
    print(f"Number of Parts (N_parts) = {n_parts}\n")

    # Step 3: Calculate the complexity.
    molecular_complexity = (n_atoms * n_bonds) / n_parts

    # Step 4: Print the final equation and result.
    print("Final Calculation:")
    print(f"Böttcher Molecular Complexity = ({n_atoms} * {n_bonds}) / {n_parts} = {int(molecular_complexity)}")

calculate_bottcher_complexity()