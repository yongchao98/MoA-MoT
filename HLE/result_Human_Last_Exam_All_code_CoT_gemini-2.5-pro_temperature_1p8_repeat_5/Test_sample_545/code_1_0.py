def get_iupac_name():
    """
    This function determines and returns the IUPAC name of the major product.

    The reaction is a thermal syn-elimination (pyrolysis) of a sulfoxide.
    Reactant: ((2-((2-methylbut-3-en-2-yl)oxy)ethyl)sulfinyl)benzene
    Structure: Ph-S(=O)-CH2(alpha)-CH2(beta)-O-C(CH3)2-CH=CH2
    The reaction eliminates Ph-SOH and forms a double bond between the alpha and beta carbons.
    Product: CH2=CH-O-C(CH3)2-CH=CH2

    IUPAC Naming (Skeletal Replacement Nomenclature):
    1. Main chain: C=C-O-C(Me)2-C=C
    2. Numbering gives lowest locant to heteroatom (O): C(1)=C(2)-O(3)-C(4)(Me)2-C(5)=C(6)
    3. Parent chain name: 3-oxahexa-1,5-diene
    4. Substituents: 4,4-dimethyl
    5. Final Name: 4,4-dimethyl-3-oxahexa-1,5-diene
    """
    
    # The final IUPAC name of the major product
    product_name = "4,4-dimethyl-3-oxahexa-1,5-diene"
    
    # Print the name as requested
    print(f"The IUPAC name of the major product is: {product_name}")

get_iupac_name()