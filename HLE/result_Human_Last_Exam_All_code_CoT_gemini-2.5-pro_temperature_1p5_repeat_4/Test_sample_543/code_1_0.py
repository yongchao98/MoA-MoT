def get_iupac_name():
    """
    This function determines and prints the IUPAC name of the reaction product.
    
    Reaction: 1,3-dibromo-2-iodobenzene + excess phenyl magnesium bromide -> ?
    
    1. The reaction is a Grignard cross-coupling where phenyl groups replace halogens.
    2. The reactivity of aryl halides is I > Br > Cl.
    3. The iodine at position 2 reacts first.
    4. Due to "excess" reagent and "reflux" conditions, the two less reactive bromines at positions 1 and 3 also react.
    5. All three halogens are substituted by phenyl groups.
    6. The final product has a central benzene ring with phenyl groups at positions 1, 2, and 3.
    7. The IUPAC name is constructed based on this structure.
    """
    
    # The parent molecule is Benzene.
    parent_molecule = "benzene"
    
    # The substituents are three phenyl groups.
    substituent_prefix = "tri"
    substituent_name = "phenyl"
    
    # The positions of the substituents on the benzene ring.
    positions = "1,2,3"
    
    # Construct the final IUPAC name.
    iupac_name = f"{positions}-{substituent_prefix}{substituent_name}{parent_molecule}"
    
    print(f"The IUPAC name of the product is: {iupac_name}")

get_iupac_name()