import sys

def get_product_name():
    """
    Determines the IUPAC name of the product from the reaction of 1,3-dibromo-2-iodobenzene
    with excess phenyl magnesium bromide.

    The reaction proceeds via substitution of all halogens by phenyl groups.
    The order of reactivity for substitution is I > Br.

    Reaction scheme:
    1. Starting material: A benzene ring with Br at pos 1, I at pos 2, and Br at pos 3.
    2. Reaction with excess PhMgBr (phenyl magnesium bromide) leads to substitution of all halogens.
       - The Iodine at position 2 is the most reactive and is substituted first.
       - The Bromines at positions 1 and 3 are substituted next.
    3. The final product is a benzene ring with phenyl groups at positions 1, 2, and 3.

    The IUPAC name for this structure is 1,2,3-triphenylbenzene.
    The numbers in the name are the locants for the phenyl groups.
    """
    
    # The parent molecule is benzene.
    parent_molecule = "benzene"
    
    # There are three phenyl substituents.
    substituent_count_prefix = "tri"
    substituent_name = "phenyl"
    
    # The positions of the substituents are 1, 2, and 3.
    locant_numbers = [1, 2, 3]
    
    # Constructing the IUPAC name
    locants_str = ",".join(map(str, locant_numbers))
    
    final_name = f"{locants_str}-{substituent_count_prefix}{substituent_name}{parent_molecule}"
    
    print(f"The IUPAC name of the product is: {final_name}")

# Execute the function to find and print the name.
get_product_name()