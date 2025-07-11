def get_product_iupac_name():
    """
    This function determines and constructs the IUPAC name of the product
    from the Oxy-Cope rearrangement shown in the image.
    """
    # 1. Define the components of the parent molecule name
    parent_ring_base = "cyclohexan"
    ketone_functional_group = "one"
    ketone_position = 1

    parent_name = f"{parent_ring_base}-{ketone_position}-{ketone_functional_group}"

    # 2. Define the components of the substituent name
    substituent_position = 3
    sub_chain_length = "but"
    sub_double_bond_position = 1
    sub_methoxy_group_position = 1
    sub_attachment_position = 3

    substituent_name = (
        f"({sub_methoxy_group_position}-methoxy"
        f"{sub_chain_length}-{sub_double_bond_position}-en"
        f"-{sub_attachment_position}-yl)"
    )

    # 3. Assemble the final IUPAC name
    final_name = f"{substituent_position}-{substituent_name}{parent_name}"

    # 4. Print the result
    print("The reaction is an Oxy-Cope rearrangement followed by keto-enol tautomerization.")
    print("The final product's IUPAC name is constructed as follows:")
    print(f"Position: {substituent_position}")
    print(f"Substituent: {substituent_name}")
    print(f"Parent Ketone: {parent_name}")
    print("\nFinal IUPAC Name:")
    print(final_name)

get_product_iupac_name()