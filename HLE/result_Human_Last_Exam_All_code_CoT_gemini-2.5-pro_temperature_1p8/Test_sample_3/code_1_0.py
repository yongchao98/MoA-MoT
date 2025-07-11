def find_largest_torsion_order():
    """
    Calculates the largest order of a non-cyclic torsion subgroup of an
    elliptic curve over the number field Q(sqrt(-3)).
    """
    
    # According to established theorems in number theory, the possible
    # non-cyclic torsion subgroups for an elliptic curve over Q(sqrt(-3)) are:
    # Z/2Z x Z/2Z, Z/2Z x Z/6Z, and Z/3Z x Z/3Z.
    # We represent these groups as tuples (m, n) for Z/mZ x Z/nZ.
    non_cyclic_groups = [
        {'m': 2, 'n': 2},
        {'m': 2, 'n': 6},
        {'m': 3, 'n': 3}
    ]

    max_order = 0
    
    print("Calculating the orders of possible non-cyclic torsion subgroups over Q(sqrt(-3)):\n")

    # Iterate through the groups, calculate their order, and find the maximum.
    for group in non_cyclic_groups:
        m = group['m']
        n = group['n']
        order = m * n
        
        # The final equation requires printing each number
        print(f"For the group Z/{m}Z x Z/{n}Z, the order is {m} * {n} = {order}")

        if order > max_order:
            max_order = order

    print("\n------------------------------------------------------------------")
    print(f"The largest possible order for such a subgroup is {max_order}.")

if __name__ == "__main__":
    find_largest_torsion_order()