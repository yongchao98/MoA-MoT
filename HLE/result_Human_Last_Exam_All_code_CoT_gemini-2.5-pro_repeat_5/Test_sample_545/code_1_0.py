def get_product_name():
    """
    This function provides the IUPAC name of the major product from the described reaction.
    
    The reaction is a tandem sequence:
    1. Sulfoxide pyrolysis (syn-elimination)
    2. Claisen rearrangement ([3,3]-sigmatropic shift)
    
    The final product is an aldehyde with the structure:
    O=CH-CH2-C(CH3)2-CH=CH2
    
    The IUPAC name is derived as:
    - Parent chain: 5 carbons with an aldehyde and a double bond -> pent-4-enal
    - Substituents: Two methyl groups at position 3 -> 3,3-dimethyl
    """
    
    # The final IUPAC name of the major product
    final_iupac_name = "3,3-dimethylpent-4-enal"
    
    # The numbers (locants) in the IUPAC name, as requested.
    # The "final equation" is interpreted as the final name.
    numbers_in_name = ['3', '3', '4']

    print(f"The IUPAC name of the major product is: {final_iupac_name}")
    
    print("\nThe numeric locants in the name are:")
    for num in numbers_in_name:
        print(num)

if __name__ == "__main__":
    get_product_name()