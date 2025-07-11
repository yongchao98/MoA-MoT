def get_product_name():
    """
    Determines and prints the IUPAC name of the product from the reaction of
    1,3-dibromo-2-iodobenzene and excess phenyl magnesium bromide.
    """
    
    # Reactants
    aryl_halide = "1,3-dibromo-2-iodobenzene"
    grignard_reagent = "excess phenyl magnesium bromide"
    
    # The reaction involves the substitution of all three halogens (I, Br, Br)
    # with phenyl groups. The order of reactivity is I > Br.
    # The final product has a central benzene ring with three phenyl substituents
    # at positions 1, 2, and 3.
    
    # IUPAC Naming
    parent_chain = "benzene"
    substituent = "phenyl"
    count = 3  # one for each halogen replaced
    prefix = "tri"
    positions = "1,2,3"
    
    product_iupac_name = f"{positions}-{prefix}{substituent}{parent_chain}"
    
    # The final equation would be C6H3BrI2 + 3 PhMgBr -> C6H3(Ph)3 + MgBrI + 2 MgBr2
    # The product is 1,2,3-triphenylbenzene. Let's output the numbers in the name as requested.
    number_1 = 1
    number_2 = 2
    number_3 = 3
    
    # Print the name ensuring the numbers are output.
    print(f"{number_1},{number_2},{number_3}-triphenylbenzene")

get_product_name()