def generate_iupac_name():
    """
    This script determines and prints the IUPAC name of the reaction product.
    The final product is 1,2,3-triphenylbenzene, which is formed by the complete
    substitution of all halogens on 1,3-dibromo-2-iodobenzene with phenyl groups
    from the excess phenyl magnesium bromide.
    """

    # Define the components of the IUPAC name based on the final product structure.
    
    # Locants (positions) of the phenyl groups on the central benzene ring.
    locants = "1,2,3"

    # Multiplicative prefix for three identical substituents.
    prefix = "tri"

    # Name of the substituent group.
    substituent = "phenyl"

    # Name of the parent ring.
    parent_molecule = "benzene"

    # Assemble the final IUPAC name according to nomenclature rules.
    # The format is locants-prefix+substituent+parent.
    final_name = f"{locants}-{prefix}{substituent}{parent_molecule}"

    print("The IUPAC name of the product is:")
    print(final_name)

generate_iupac_name()