def solve_chemistry_problem():
    """
    This function determines the IUPAC name of the product from the reaction scheme
    and prints it. The name is constructed from its components.
    """

    # Components of the IUPAC name
    ester_alkyl_group = "ethyl"
    ring_system_parent = "thiophene"
    saturation_locants = [2, 5]
    saturation_prefix = "dihydro"
    substituent_locant = 3
    substituent_suffix = "carboxylate"

    # Assemble the final IUPAC name
    product_name = (
        f"{ester_alkyl_group} "
        f"{saturation_locants[0]},{saturation_locants[1]}-"
        f"{saturation_prefix}{ring_system_parent}-"
        f"{substituent_locant}-{substituent_suffix}"
    )

    print(product_name)

solve_chemistry_problem()