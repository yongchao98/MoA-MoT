import re

def solve_wittig_reaction():
    """
    Solves for the product of a Wittig reaction between pivalaldehyde
    and (2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane.
    """

    # 1. Define the reactants and their key fragments for the reaction.
    aldehyde = {
        "name": "pivalaldehyde",
        "r_group_fragment": "(tert-butyl)" # Represents (CH3)3C-
    }

    ylide = {
        "name": "(2-(2-chlorophenyl)ethylidene)triphenyl-l5-phosphane",
        "alkylidene_fragment": "CH-CH2-(2-chlorophenyl)", # The part that replaces the oxygen
        "type": "non-stabilized"  # Determines stereochemistry
    }

    # 2. Determine the product's stereochemistry based on the ylide type.
    # Non-stabilized ylides favor the (Z)-isomer.
    if ylide["type"] == "non-stabilized":
        stereochemistry = "(Z)-"
    else:
        stereochemistry = "(E)-"

    # 3. Assemble the IUPAC name for the product.
    # This process mimics the manual IUPAC naming rules.
    # Structure: (CH3)3C-CH=CH-CH2-(2-chlorophenyl)
    # Longest chain containing C=C is 5 carbons long (pentene).
    # Double bond is at position 2 (pent-2-ene).
    # Substituents are 1-(2-chlorophenyl) and 4,4-dimethyl.
    
    parent_chain = "pent"
    double_bond_locant = "2"
    substituent1 = "1-(2-chlorophenyl)"
    substituent2 = "4,4-dimethyl"

    # Combine parts to form the full name, which we will treat as the "final equation"
    product_name = f"{stereochemistry}{substituent1}-{substituent2}{parent_chain}-{double_bond_locant}-ene"
    
    print(f"The major product of the Wittig reaction is:")
    print(product_name)
    print()

    # 4. As requested, output each number from the "final equation" (the IUPAC name).
    # We use regex to find all digits in the product name string.
    numbers_in_name = re.findall(r'\d', product_name)
    
    print("The numbers in the final product name are:")
    # Print each number found in the name.
    for number in numbers_in_name:
        print(number)

solve_wittig_reaction()