import re

def wittig_reaction_product():
    """
    This function outlines the Wittig reaction between pivalaldehyde and a specific ylide,
    determines the product, and prints the details as requested.
    """

    # --- Reactants and Products ---
    reactant_aldehyde = "Pivalaldehyde [(CH3)3C-CHO]"
    reactant_ylide = "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane [(2-Cl-C6H4)-CH2-CH=P(C6H5)3]"
    product_alkene_name = "(Z)-1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
    byproduct = "Triphenylphosphine oxide [O=P(C6H5)3]"

    # --- Print the Reaction Equation ---
    print("The Wittig reaction between the specified reactants is as follows:")
    print(f"Reactant 1: {reactant_aldehyde}")
    print(f"Reactant 2: {reactant_ylide}")
    print("\n--->\n")
    print(f"Major Product: {product_alkene_name}")
    print(f"Byproduct: {byproduct}")
    print("-" * 50)

    # --- Output each number in the final equation (product name) ---
    print("\nAs requested, here are the numbers present in the IUPAC name of the final product:")

    # Find all sequences of digits in the product's name
    numbers_in_name = re.findall(r'\d+', product_alkene_name)

    # Print each number found
    print(f"The numbers are: {', '.join(numbers_in_name)}")

    # Explain the meaning of each number
    print("\nWhere these numbers signify:")
    print("  - 1: The position of the '(2-chlorophenyl)' substituent on the main pentene chain.")
    print("  - 2: The position of the 'chloro' substituent on the phenyl ring.")
    print("  - 4, 4: The positions of the two 'dimethyl' substituents on the pentene chain.")
    print("  - 2: The starting position of the 'ene' (double bond) on the pentene chain.")


# Execute the function
wittig_reaction_product()
