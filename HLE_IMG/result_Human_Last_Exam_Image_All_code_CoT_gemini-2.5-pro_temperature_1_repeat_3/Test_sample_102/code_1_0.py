def get_iupac_name():
    """
    This function determines and prints the IUPAC name of the reaction product.
    """
    # 1. Identify the parent chain based on the principal functional group (aldehyde)
    parent_chain_length = 3
    parent_name = "propanal"
    print(f"Step 1: The principal functional group is an aldehyde.")
    print(f"Step 2: The longest carbon chain containing the aldehyde is {parent_chain_length} carbons long.")
    print(f"Step 3: The parent name is '{parent_name}'.\n")

    # 2. Identify the substituents and their locations
    substituents = {
        "methoxy": 2,
        "methyl": 3,
        "cyclohex-2-en-1-yl": 3
    }
    
    print("Step 4: Identify the substituents on the parent chain:")
    for sub, loc in substituents.items():
        print(f"  - A '{sub}' group is at position {loc}.")
    print("")

    # 3. Assemble the full IUPAC name (substituents are ordered alphabetically)
    # Alphabetical order: cyclohexenyl, methoxy, methyl
    final_name = f"{substituents['cyclohex-2-en-1-yl']}-({substituents['cyclohex-2-en-1-yl']})-{substituents['methoxy']}-{substituents['methoxy']}-{substituents['methyl']}-{substituents['methyl']}{parent_name}"
    # Corrected assembly
    final_name = f"3-(cyclohex-2-en-1-yl)-2-methoxy-3-methyl{parent_name}"
    
    print("Step 5: Assemble the final name, ordering substituents alphabetically.")
    print("Final IUPAC Name: " + final_name)
    print("\nThe numbers used in the final name are: 3, 2, 1, 2, 3.")


get_iupac_name()