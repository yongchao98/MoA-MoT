def get_iupac_name():
    """
    This function constructs and prints the IUPAC name of the final product
    based on the chemical analysis.
    """

    # The final product's structure is (CH3)2C=CH-CH2-CH2-CHO.
    # We will derive its IUPAC name.

    prefix_number = '5'
    prefix_number_2 = '5'
    prefix_substituent = 'dimethyl'
    parent_chain = 'pent'
    double_bond_locant = '4'
    suffix_double_bond = 'en'
    suffix_functional_group = 'al'

    # The IUPAC name is constructed by combining these parts.
    # The rule is: (locants)-(substituents)(parent chain)-(locant)-(suffix)

    final_name = f"{prefix_number},{prefix_number_2}-{prefix_substituent}{parent_chain}-{double_bond_locant}-{suffix_functional_group}"

    print("The IUPAC name for the major product is derived as follows:")
    print(f"Substituent locants (positions): {prefix_number}, {prefix_number_2}")
    print(f"Substituent name: {prefix_substituent}")
    print(f"Parent chain length: {parent_chain} (5 carbons)")
    print(f"Double bond locant: {double_bond_locant}")
    print(f"Functional group suffix: {suffix_double_bond}{suffix_functional_group}")

    print("\n" + "="*40)
    print("Final IUPAC Name:")
    print(final_name)
    print("="*40)

if __name__ == "__main__":
    get_iupac_name()