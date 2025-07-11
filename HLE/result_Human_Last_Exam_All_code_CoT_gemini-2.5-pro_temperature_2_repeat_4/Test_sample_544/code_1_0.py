import sys

def get_product_name():
    """
    This function determines the IUPAC name of the product from the reaction of
    methyl phenyl sulfoxide with triflic anhydride and trimethylsilyl cyanide.
    
    The reaction is a Pummerer rearrangement.
    1.  Reactants: methyl phenyl sulfoxide, triflic anhydride, trimethylsilyl cyanide.
    2.  The sulfoxide is activated by the anhydride.
    3.  A Pummerer rearrangement occurs, forming a thionium ion intermediate: [Ph-S=CH2]+.
    4.  The cyanide nucleophile attacks the electrophilic carbon.
    5.  The final product structure is Ph-S-CH2-CN.
    6.  The IUPAC name for Ph-S-CH2-CN is 2-(phenylthio)acetonitrile.
    """
    
    # The parent chain is acetonitrile (-CH2-CN)
    parent_chain = "acetonitrile"
    
    # The substituent is a phenylthio group (Ph-S-)
    substituent = "phenylthio"
    
    # The position of the substituent on the acetonitrile chain is 2
    position = 2
    
    # The final name is constructed according to IUPAC rules.
    final_name = f"{position}-({substituent}){parent_chain}"
    
    return final_name

if __name__ == "__main__":
    product_name = get_product_name()
    print("The IUPAC name of the reaction product is:")
    print(product_name)
