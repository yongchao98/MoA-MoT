def generate_iupac_name():
    """
    This function generates the IUPAC name for the final reaction product
    based on its determined chemical structure: OHC-CH2-CH2-CH=C(CH3)2
    """

    # IUPAC Naming Rule Application
    # 1. Principal Functional Group: Aldehyde -> suffix 'al'. Aldehyde carbon is C1.
    # 2. Longest Carbon Chain (Parent): 5 carbons including the aldehyde -> root 'pent'.
    # 3. Unsaturation (Double Bond): Starts at C4 -> infix '-4-en'.
    # 4. Substituents (Prefix): Two methyl groups on C5 -> prefix '5,5-dimethyl'.

    # Component parts of the name
    prefix_locants = "5,5"
    prefix_text = "dimethyl"
    root = "pent"
    unsaturation_locant = 4
    unsaturation_type = "en"
    suffix = "al"

    # Assemble the final name according to IUPAC preferred nomenclature
    final_name = f"{prefix_locants}-{prefix_text}{root}-{unsaturation_locant}-{unsaturation_type}{suffix}"

    print("The final product's IUPAC name is determined from its structure:")
    print("Structure: OHC-C(2)H2-C(3)H2-C(4)H=C(5)(CH3)2")
    print("\nThe numerical locants in the name are:")
    print(f"Methyl group locants: {prefix_locants.split(',')[0]} and {prefix_locants.split(',')[1]}")
    print(f"Double bond locant: {unsaturation_locant}")

    print("\nFinal Assembled IUPAC Name:")
    print(final_name)

generate_iupac_name()