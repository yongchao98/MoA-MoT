def solve_chemistry_problem():
    """
    This script determines the IUPAC name of the major product from the described reaction.
    The reaction involves a tandem sulfoxide elimination followed by a Claisen rearrangement.
    """

    # Step 1: Define the structure of the final product based on reaction analysis.
    # Reactant: Ph-S(=O)-CH2-CH2-O-C(CH3)2-CH=CH2
    # Intermediate (via sulfoxide elimination): CH2=CH-O-C(CH3)2-CH=CH2
    # Final Product (via Claisen rearrangement): O=CH-CH2-CH2-CH=C(CH3)2
    product_structure = "O=CH-CH2-CH2-CH=C(CH3)2"

    # Step 2: Systematically determine the IUPAC name for the product.
    
    # Substituents: Two methyl groups on carbon 5.
    substituents = "5,5-dimethyl"
    
    # Parent chain length: 6 carbons including the aldehyde.
    parent_chain = "hex"
    
    # Position of the double bond: Between carbons 4 and 5.
    unsaturation_locator = "4"
    unsaturation_type = "en"
    
    # Principal functional group: Aldehyde.
    functional_group = "al"
    
    # Step 3: Assemble the final IUPAC name.
    # The format is: [substituents][parent_chain]-[locator]-[unsaturation][suffix]
    final_iupac_name = f"{substituents}{parent_chain}-{unsaturation_locator}-{unsaturation_type}{functional_group}"

    # Print the result.
    print("The IUPAC name of the major product is:")
    print(final_iupac_name)

if __name__ == "__main__":
    solve_chemistry_problem()