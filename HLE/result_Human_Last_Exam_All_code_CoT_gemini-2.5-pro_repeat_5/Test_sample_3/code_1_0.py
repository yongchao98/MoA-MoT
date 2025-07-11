def find_largest_torsion_order():
    """
    Calculates the largest order of a non-cyclic torsion subgroup of an
    elliptic curve over the number field Q(sqrt(-3)).

    The calculation is based on the known classification of possible
    torsion structures over this field.
    """
    # A list of tuples (m, n) representing the non-cyclic groups Z/mZ x Z/nZ.
    # The structure of these groups is a known result from number theory.
    non_cyclic_groups = []

    # Add groups of the form Z/2Z x Z/2NZ for N in {1, ..., 7}
    for N in range(1, 8):
        non_cyclic_groups.append((2, 2 * N))

    # Add groups of the form Z/3Z x Z/3NZ for N in {1, 2, 3}
    for N in range(1, 4):
        non_cyclic_groups.append((3, 3 * N))

    # Add the group Z/4Z x Z/4Z
    non_cyclic_groups.append((4, 4))

    # Initialize variables to store the maximum order and the corresponding group
    max_order = 0
    
    print("Calculating the order for each possible non-cyclic torsion subgroup:")
    print("-----------------------------------------------------------------")
    
    # Iterate through the list of groups, calculate their orders, and find the max
    for m, n in non_cyclic_groups:
        order = m * n
        print(f"Group Z/{m}Z x Z/{n}Z has order {m} * {n} = {order}")
        if order > max_order:
            max_order = order

    print("\n-----------------------------------------------------------------")
    print(f"The largest order among these non-cyclic torsion subgroups is {max_order}.")

if __name__ == "__main__":
    find_largest_torsion_order()