def get_diels_alder_product():
    """
    This function determines the product of the reaction between
    butadiene and 1,1-dichloro-2,2-difluoroethene and prints the result.
    """
    
    # Define the reactants
    diene = "Butadiene"
    dienophile = "1,1-dichloro-2,2-difluoroethene"
    
    # The reaction is a Diels-Alder [4+2] cycloaddition, forming a cyclohexene ring.
    
    # Determine the substituents and their locants (positions) on the final product ring.
    # IUPAC rules state we number the ring to give the double bond positions 1 and 2.
    # Then, we assign the lowest possible numbers to the substituents.
    # When locant sets are tied (both are 4,4,5,5), the alphabetically first
    # substituent ('chloro') gets the lower number.
    
    locant_for_chloro = "4,4"
    substituent_chloro = "dichloro"
    
    locant_for_fluoro = "5,5"
    substituent_fluoro = "difluoro"
    
    parent_ring = "cyclohexene"

    # Assemble the final product name
    product_name = f"{locant_for_chloro}-{substituent_chloro}-{locant_for_fluoro}-{substituent_fluoro}{parent_ring}"

    # Print the explanation and the final result
    print(f"The reaction between {diene} and {dienophile} is a Diels-Alder reaction.")
    print("This [4+2] cycloaddition forms a six-membered ring.")
    print("\nThe final product is named based on its structure:")
    print(f" - The numbers {locant_for_chloro} are the positions for the two chlorine atoms.")
    print(f" - The numbers {locant_for_fluoro} are the positions for the two fluorine atoms.")
    print("\nFinal Product Name:")
    print(product_name)

get_diels_alder_product()