def get_reaction_product():
    """
    This function describes the product of the Diels-Alder reaction between
    butadiene and 1,1-dichloro-2,2-difluoroethene.
    """
    
    # Define the reactants and the product by their names
    diene = "1,3-Butadiene"
    dienophile = "1,1-dichloro-2,2-difluoroethene"
    product_name = "4,4-dichloro-5,5-difluorocyclohexene"
    product_formula = "C6H6Cl2F2"

    # The reaction is a Diels-Alder [4+2] cycloaddition
    reaction_type = "Diels-Alder [4+2] cycloaddition"

    # Print the reaction details
    print(f"The reaction between {diene} and {dienophile} is a {reaction_type}.")
    print("\nThe overall reaction equation is:")
    print(f"{diene} + {dienophile} -> {product_name}")
    
    print("\n--- Product Details ---")
    print(f"IUPAC Name: {product_name}")
    print(f"Molecular Formula: {product_formula}")

    # As per the instruction to output each number in the final equation,
    # we highlight the locant numbers used in the product's name.
    locant_numbers = "4, 4, 5, 5"
    print(f"The locant numbers in the product's name are: {locant_numbers}.")

# Execute the function to print the result
get_reaction_product()