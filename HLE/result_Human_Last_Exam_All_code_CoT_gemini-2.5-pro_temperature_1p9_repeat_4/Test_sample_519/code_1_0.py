def solve_properties():
    """
    This function calculates the properties of the three given categories
    and prints them in the specified format.
    """

    # Properties for X1
    dim1 = 3 * 11
    profile1 = f"[S, s, {dim1}]"

    # Properties for X2
    dim2 = 4 - 1
    profile2 = f"[DM, s, irr, {dim2}]"

    # Properties for X3
    genus_g = 7
    dim3 = genus_g + 1
    profile3 = f"[A, s, uc, {dim3}]"

    # Combine and print the results
    final_answer = f"{profile1} {profile2} {profile3}"
    print(final_answer)

solve_properties()