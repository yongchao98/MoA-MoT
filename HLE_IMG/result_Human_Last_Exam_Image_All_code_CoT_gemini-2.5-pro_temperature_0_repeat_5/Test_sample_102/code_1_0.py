def get_product_name():
    """
    This function determines and prints the IUPAC name of the product
    from the given Oxy-Cope rearrangement reaction.
    """
    # Locants (position numbers) in the final IUPAC name
    locant_substituent_on_ring = 6
    locant_methoxy_on_side_chain = 1
    locant_double_bond_on_side_chain = 1
    locant_attachment_point_on_side_chain = 3

    # Stereochemical descriptor for the side chain double bond
    alkene_geometry = "Z"

    # Parts of the name
    parent_ring = "cyclohexanone"
    side_chain_base = "methoxybut"
    side_chain_ene = "en"
    side_chain_yl = "yl"

    # Constructing the substituent name
    substituent_name = f"(({alkene_geometry})-{locant_methoxy_on_side_chain}-{side_chain_base}-{locant_double_bond_on_side_chain}-{side_chain_ene}-{locant_attachment_point_on_side_chain}-{side_chain_yl})"

    # Constructing the final IUPAC name
    final_name = f"{locant_substituent_on_ring}-{substituent_name}{parent_ring}"

    # Print the numbers used in the name as per the instructions
    print(f"The reaction is an Oxy-Cope rearrangement followed by keto-enol tautomerization.")
    print(f"The final product is a substituted cyclohexanone.")
    print(f"Numbering for the final IUPAC name:")
    print(f"Position of the substituent on the cyclohexanone ring: {locant_substituent_on_ring}")
    print(f"Position of the methoxy group in the side chain: {locant_methoxy_on_side_chain}")
    print(f"Position of the double bond in the side chain: {locant_double_bond_on_side_chain}")
    print(f"Attachment point of the side chain (as a radical): {locant_attachment_point_on_side_chain}")
    print("-" * 20)
    print("Final IUPAC Name:")
    print(final_name)

get_product_name()