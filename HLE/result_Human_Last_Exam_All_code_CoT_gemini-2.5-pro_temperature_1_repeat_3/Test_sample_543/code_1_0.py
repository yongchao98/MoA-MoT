def solve_reaction():
    """
    This function determines and prints the IUPAC name of the product from the reaction of
    1,3-dibromo-2-iodobenzene with excess phenyl magnesium bromide.
    """
    # The reaction involves the substitution of all three halogen atoms (I, Br, Br)
    # on the benzene ring with phenyl groups from the excess Grignard reagent (phenyl magnesium bromide).
    # The order of reactivity is I > Br > Cl.
    # 1. Iodine at position 2 is replaced by a phenyl group.
    # 2. Bromine at position 1 is replaced by a phenyl group.
    # 3. Bromine at position 3 is replaced by a phenyl group.
    # The final product is a benzene ring with three phenyl groups at positions 1, 2, and 3.

    product_name = "1,2,3-triphenylbenzene"
    
    # Printing the final IUPAC name
    print(f"The IUPAC name of the final product is: {product_name}")

solve_reaction()