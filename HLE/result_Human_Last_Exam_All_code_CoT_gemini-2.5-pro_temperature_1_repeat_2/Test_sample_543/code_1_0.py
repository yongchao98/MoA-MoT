def get_product_iupac_name():
    """
    Determines and prints the IUPAC name of the product from the reaction of
    1,3-dibromo-2-iodobenzene and excess phenyl magnesium bromide.
    """
    # 1. Define the reactants and reaction conditions
    aryl_halide = "1,3-dibromo-2-iodobenzene"
    grignard_reagent = "phenyl magnesium bromide"
    condition = "excess"

    print(f"Reaction: {aryl_halide} + {condition} {grignard_reagent}")
    print("-" * 50)

    # 2. Explain the reaction logic
    print("This is a Grignard cross-coupling reaction. The phenyl group from the Grignard reagent")
    print("replaces the halogen atoms on the benzene ring.")
    print("\nThe reactivity order of halogens as leaving groups is I > Br > Cl.")
    print(f"Therefore, the iodine at position 2 of {aryl_halide} reacts first,")
    print("followed by the two bromine atoms at positions 1 and 3.")
    print("\nSince the Grignard reagent is in excess, all three halogens will be replaced by phenyl groups.")
    print("-" * 50)

    # 3. Determine and construct the IUPAC name of the final product
    print("The final product is a benzene ring substituted with three phenyl groups.")
    
    # Define components of the IUPAC name
    locants = "1,2,3"
    prefix = "tri"
    substituent = "phenyl"
    parent_chain = "benzene"

    # Assemble the final name
    final_product_name = f"{locants}-{prefix}{substituent}{parent_chain}"

    print("\nConstructing the IUPAC name:")
    print(f"  - Parent molecule: {parent_chain}")
    print(f"  - Substituents: Three '{substituent}' groups, indicated by the prefix '{prefix}'.")
    print(f"  - Locants (positions): The numbers {locants.replace(',', ', ')} indicate the attachment points.")
    
    print("\nFinal IUPAC Name of the Product:")
    # The final output prints the name including each number as requested.
    print(final_product_name)

# Run the function to get the answer
get_product_iupac_name()