import math

def calculate_bottcher_complexity():
    """
    Calculates the Böttcher Molecular Complexity for the product of the
    Favorskii rearrangement of 2-chlorocyclohexanone.
    """

    # Step 1: Identify the product and its formula.
    # The Favorskii rearrangement of 2-chlorocyclohexanone yields cyclopentanecarboxylic acid.
    # The molecular formula is C6H10O2.
    product_name = "Cyclopentanecarboxylic acid"
    product_formula = "C6H10O2"

    print(f"The product of the reaction is {product_name} ({product_formula}).")
    print("-" * 30)

    # Step 2: Determine the molecular properties.
    # Number of atoms: 6 (Carbon) + 10 (Hydrogen) + 2 (Oxygen)
    N_atoms = 6 + 10 + 2

    # Number of bonds (connections):
    # C-C in ring: 5
    # C-C to carboxyl: 1
    # C-H on ring: 9
    # C=O bond: 1
    # C-O bond: 1
    # O-H bond: 1
    # Total = 5 + 1 + 9 + 1 + 1 + 1 = 18
    N_bonds = 18

    # Number of chiral centers: The carbon attached to the -COOH group is not
    # chiral because the paths around the ring are identical.
    N_chiral_centers = 0

    # Number of parts: The molecule is a single connected entity.
    N_parts = 1

    print("Calculating the necessary parameters:")
    print(f"Number of atoms (N_atoms) = {N_atoms}")
    print(f"Number of bonds (N_bonds) = {N_bonds}")
    print(f"Number of chiral centers (N_chiral_centers) = {N_chiral_centers}")
    print(f"Number of parts (N_parts) = {N_parts}")
    print("-" * 30)

    # Step 3: Define and apply the complexity formula.
    # Formula: (N_atoms * N_bonds * (N_chiral_centers + 1)) / N_parts
    print("Using the formula: Complexity = (N_atoms * N_bonds * (N_chiral_centers + 1)) / N_parts")
    
    # Calculate the complexity
    complexity = (N_atoms * N_bonds * (N_chiral_centers + 1)) / N_parts
    
    # Step 4: Display the final equation and result.
    print("\nFinal Calculation:")
    print(f"Complexity = ({N_atoms} * {N_bonds} * ({N_chiral_centers} + 1)) / {N_parts}")
    print(f"Complexity = {int(complexity)}")


if __name__ == "__main__":
    calculate_bottcher_complexity()