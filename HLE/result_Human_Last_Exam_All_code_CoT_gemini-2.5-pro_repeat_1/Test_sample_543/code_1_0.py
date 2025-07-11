def get_iupac_name():
    """
    This function determines and prints the IUPAC name of the reaction product.
    
    The reaction involves 1,3-dibromo-2-iodobenzene reacting with excess
    phenyl magnesium bromide. The reactivity of halogens in this Grignard
    cross-coupling reaction is I > Br. The excess reagent ensures all
    halogens are substituted.
    
    Step 1: The iodine at position 2 is replaced by a phenyl group.
    Step 2: The bromines at positions 1 and 3 are replaced by phenyl groups.
    
    The final product has a central benzene ring with three phenyl groups
    attached at positions 1, 2, and 3.
    """
    # The numbers for the locants in the final IUPAC name
    locant_1 = 1
    locant_2 = 2
    locant_3 = 3
    
    # The substituent name and parent chain
    prefix = "tri"
    substituent = "phenyl"
    parent_chain = "benzene"
    
    # Constructing the final IUPAC name string
    # The problem asks to output each number in the final name,
    # so we print them as part of the full name.
    final_name = f"{locant_1},{locant_2},{locant_3}-{prefix}{substituent}{parent_chain}"
    
    print(final_name)

get_iupac_name()