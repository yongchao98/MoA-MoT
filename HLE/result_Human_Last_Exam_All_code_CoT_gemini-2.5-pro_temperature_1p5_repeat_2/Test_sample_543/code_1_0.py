def get_product_name():
    """
    This function determines and prints the IUPAC name of the product from the described reaction.

    Reaction: 1,3-dibromo-2-iodobenzene reacts with excess phenyl magnesium bromide
              in refluxing THF, followed by an aqueous work-up.

    Analysis:
    1.  The reactivity of the halogens is I > Br. The most reactive C-I bond is substituted first.
    2.  Excess Grignard reagent and reflux conditions ensure the less reactive C-Br bonds also react.
    3.  All three halogens are replaced by phenyl groups.
    4.  The final product has a benzene core with phenyl groups at positions 1, 2, and 3.
    """
    
    # Define the numbers for the positions of the phenyl groups
    position_1 = 1
    position_2 = 2
    position_3 = 3
    
    # Construct the IUPAC name
    product_name = f"{position_1},{position_2},{position_3}-triphenylbenzene"
    
    print(f"The IUPAC name of the final product is: {product_name}")

# Execute the function to print the name
get_product_name()