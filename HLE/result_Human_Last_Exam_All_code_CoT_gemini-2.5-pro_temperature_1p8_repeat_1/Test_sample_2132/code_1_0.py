import math

def calculate_bottcher_complexity():
    """
    Calculates the Böttcher Molecular Complexity for the product of the
    Favorskii rearrangement of 2-chlorocyclohexanone.
    """

    # Step 1: Identify the product and its formula.
    # The reaction is the Favorskii rearrangement of 2-chlorocyclohexanone.
    # The product is cyclopentanecarboxylic acid.
    # The chemical formula for cyclopentanecarboxylic acid is C6H10O2.
    product_name = "cyclopentanecarboxylic acid"
    product_formula = "C6H10O2"
    num_carbons = 6
    num_hydrogens = 10
    num_oxygens = 2

    print("This script calculates the Böttcher Molecular Complexity for the product of the Favorskii rearrangement of 2-chlorocyclohexanone.")
    print(f"The product of the reaction is {product_name}, with the chemical formula {product_formula}.\n")
    print("The Böttcher Molecular Complexity is calculated using the formula: MC = (N_atoms * N_bonds) / N_parts\n")

    # Step 2: Calculate the parameters N_atoms, N_bonds, and N_parts.

    # N_atoms is the total number of atoms.
    N_atoms = num_carbons + num_hydrogens + num_oxygens
    print(f"Calculating N_atoms (total number of atoms):")
    print(f"N_atoms = {num_carbons} (C) + {num_hydrogens} (H) + {num_oxygens} (O) = {N_atoms}\n")

    # N_bonds is the total number of valence bonds (double bonds count as 2).
    # This is calculated as (sum of valences) / 2.
    # Valence: C=4, H=1, O=2.
    sum_of_valences = (num_carbons * 4) + (num_hydrogens * 1) + (num_oxygens * 2)
    N_bonds = sum_of_valences // 2
    print(f"Calculating N_bonds (total number of covalent bonds):")
    print(f"Sum of valences = ({num_carbons} * 4) + ({num_hydrogens} * 1) + ({num_oxygens} * 2) = {sum_of_valences}")
    print(f"N_bonds = {sum_of_valences} / 2 = {N_bonds}\n")

    # N_parts is the number of cyclic fragments.
    # Cyclopentanecarboxylic acid has one ring.
    N_parts = 1
    print(f"Determining N_parts (number of rings):")
    print(f"The molecule has one cyclopentane ring, so N_parts = {N_parts}\n")

    # Step 3: Calculate the molecular complexity.
    molecular_complexity = (N_atoms * N_bonds) / N_parts

    # Step 4: Output the final calculation and result.
    print("Final Calculation:")
    print(f"Böttcher Molecular Complexity = ({N_atoms} * {N_bonds}) / {N_parts}")
    print(f"Result = {molecular_complexity}")


if __name__ == '__main__':
    calculate_bottcher_complexity()
<<<342.0>>>