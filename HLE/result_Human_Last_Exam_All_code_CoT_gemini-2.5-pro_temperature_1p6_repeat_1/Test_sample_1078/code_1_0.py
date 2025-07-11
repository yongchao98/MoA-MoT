def count_weyl_groups_rank_4():
    """
    This function determines the number of non-isomorphic finite Weyl groups of rank 4.

    Finite Weyl groups are classified by irreducible root systems, which correspond to
    the connected Dynkin diagrams: A_n, B_n, C_n, D_n, E_6, E_7, E_8, F_4, G_2.
    The rank of the group is the number of nodes in its diagram. We need to find
    all such diagrams with exactly 4 nodes.
    """

    # Dictionary to hold the count for each type.
    # We count 1 for each type that has a rank-4 representative.
    # Note: B_n and C_n have isomorphic Weyl groups, so they are considered as one type.
    contributions = {
        "A_4": 1,        # The group from the A_n series with n=4
        "B_4 / C_4": 1,  # The group from B_n/C_n series with n=4 (W(B_4) is isomorphic to W(C_4))
        "D_4": 1,        # The group from the D_n series with n=4
        "F_4": 1,        # The exceptional group F_4 is of rank 4 by definition
    }

    print("The non-isomorphic finite Weyl groups of rank 4 are derived from the following Dynkin diagrams:")
    for group_type, count in contributions.items():
        print(f"- Type {group_type}: {count} group")

    # Calculate the total by summing the contributions
    total_count = sum(contributions.values())
    
    # Get the individual numbers for the final equation
    numbers = list(contributions.values())
    
    print("\nTo find the total number, we sum the counts from each type:")
    # Creates a string like "1 + 1 + 1 + 1 = 4"
    equation_str = " + ".join(map(str, numbers)) + f" = {total_count}"
    
    print(equation_str)


# Execute the function to print the result
count_weyl_groups_rank_4()