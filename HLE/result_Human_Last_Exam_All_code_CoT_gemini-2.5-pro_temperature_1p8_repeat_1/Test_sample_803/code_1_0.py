def find_nonabelian_filled_groups():
    """
    Identifies and describes the nonabelian filled groups of order 2*q^m,
    based on a specific characterization from group theory literature.

    Note: The definition of a "filled group" has led to conflicting
    characterizations in research. This solution follows the characterization that
    a group is filled if and only if all of its subgroups have cyclic Sylow
    subgroups. This leads to a constructive result.
    """

    print("--- General Solution for nonabelian filled groups of order 2*q^m ---")
    print("For an odd prime 'q' and natural number 'm', the only family of nonabelian groups of order 2*q^m that are 'filled' (under the chosen definition) are the Dihedral Groups.\n")

    print("1. Identification of the Group:")
    print("The group is the Dihedral Group of order 2*q^m, denoted D_{2q^m}.\n")

    print("2. Group Presentation:")
    print("The standard presentation for this group is given by generators 'r' (for rotation) and 's' (for reflection) and a set of relations:")
    print("G = < r, s | r^(q^m) = 1, s^2 = 1, srs = r^-1 >\n")

    print("3. Numbers in the Final Equations:")
    print("As requested, here are the numbers that appear in the defining equations of the presentation:")
    print("- In 'r^(q^m) = 1': The exponent of 'r' is the number q^m. The identity element '1' is also part of this equation.")
    print("- In 's^2 = 1': The exponent of 's' is the number 2.")
    print("- In 'srs = r^-1': The exponent of 'r' on the right side is the number -1. The exponents of r and s on the left are implicitly 1.\n")

    print("--- Example: q = 3, m = 2 ---")
    q = 3
    m = 2
    order = 2 * (q**m)
    r_exponent = q**m

    print(f"For q={q} and m={m}, the order of the group is 2 * {q}^{m} = {order}.")
    print(f"The group is the Dihedral Group D_{{{order}}}.\n")

    print("Its presentation is:")
    print(f"G = < r, s | r^{r_exponent} = 1, s^2 = 1, srs = r^-1 >\n")

    print("The numbers in these specific equations are:")
    print(f"- From 'r^{r_exponent} = 1': The number is {r_exponent}.")
    print("- From 's^2 = 1': The number is 2.")
    print("- From 'srs = r^-1': The number is -1.")

# Execute the function to display the solution.
find_nonabelian_filled_groups()
