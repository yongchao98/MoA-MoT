def solve_scheme_properties():
    """
    This function prints the lexicographically sorted list of all maximal subsets of
    properties {A,B,C,D,E} for which a scheme X exists.
    The derivation of these subsets is based on principles of algebraic geometry.
    """
    # The maximal subsets, determined through mathematical analysis and sorted lexicographically.
    maximal_subsets = [
        "{A,B}",
        "{A,C,D}",
        "{A,C,E}",
        "{B,D}"
    ]

    # Join the formatted sets into a single string as per the example output format.
    final_answer = ", ".join(maximal_subsets)

    print(final_answer)

solve_scheme_properties()