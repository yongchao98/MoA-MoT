def find_largest_torsion_order():
    """
    This function determines the largest order of a non-cyclic torsion subgroup
    of an elliptic curve over the number field Q(sqrt(-3)).

    The calculation is based on the complete classification of possible torsion
    subgroups for elliptic curves over quadratic fields.
    """

    # According to the classification, the possible non-cyclic torsion subgroups
    # G = Z/mZ x Z/nZ over Q(sqrt(-3)) are known. We list them below,
    # represented by the pair (m, n).
    possible_groups = []

    # Groups of the form Z/2Z x Z/2NZ for N in {1, 2, 3, 4, 5, 6}
    for N in range(1, 7):
        possible_groups.append((2, 2 * N))

    # Groups of the form Z/3Z x Z/3NZ for N in {1, 2}
    for N in range(1, 3):
        possible_groups.append((3, 3 * N))

    # Group of the form Z/4Z x Z/4Z
    possible_groups.append((4, 4))

    # Group of the form Z/6Z x Z/6Z
    possible_groups.append((6, 6))

    max_order = 0
    group_with_max_order = (0, 0)

    # Iterate through the list to find the group with the maximum order
    for m, n in possible_groups:
        order = m * n
        if order > max_order:
            max_order = order
            group_with_max_order = (m, n)

    m, n = group_with_max_order
    print(f"The non-cyclic torsion subgroup with the largest order is Z/{m}Z x Z/{n}Z.")
    print("The calculation for its order is:")
    # Final output showing the numbers in the equation
    print(f"{m} * {n} = {max_order}")
    print(f"\nTherefore, the largest order of a non-cyclic torsion subgroup is {max_order}.")

if __name__ == '__main__':
    find_largest_torsion_order()