import re

def solve_wittig_reaction():
    """
    This function outlines the Wittig reaction and prints the final product's name.
    """
    # Define the reactants and products by name for the final equation
    aldehyde = "pivalaldehyde"
    ylide = "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane"
    # The major product is the Z-isomer
    product_alkene = "(Z)-1-(2-chlorophenyl)-4,4-dimethylpent-2-ene"
    byproduct = "triphenylphosphine oxide"

    # The full reaction equation
    reaction_equation = f"{aldehyde} + {ylide} -> {product_alkene} + {byproduct}"

    print("The final reaction equation is:")
    print(reaction_equation)
    print("\n--------------------------------------------------")
    print("The name of the main organic product is:")
    print(product_alkene)

    # As requested, output each number present in the final product's IUPAC name.
    # We use regex to find all numbers in the product name string.
    numbers_in_name = re.findall(r'\d', product_alkene)
    
    print("\nThe numbers in the final product's IUPAC name are:")
    # The numbers specify the positions of substituents and the double bond.
    # 1: position of the (2-chlorophenyl) group
    # 2: position of the chloro group on the phenyl ring
    # 4,4: positions of the two methyl groups
    # 2: position of the pentene double bond
    print(', '.join(numbers_in_name))

solve_wittig_reaction()