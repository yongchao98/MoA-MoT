import sys

def solve_torsion_order():
    """
    Calculates the largest order of a non-cyclic torsion subgroup of an elliptic curve
    over the number field Q(sqrt(-3)).
    
    The possible torsion structures are known from deep theorems in number theory.
    This script uses the known classification for Q(sqrt(-3)).
    """

    # According to the classification of torsion subgroups of elliptic curves over
    # quadratic fields, the non-cyclic torsion subgroups that can occur over K = Q(sqrt(-3))
    # are of the form Z/mZ x Z/nZ, where (m, n) are in the following list.
    possible_groups = [
        # Groups of the form Z/2Z x Z/2NZ for N in {1, 2, 3, 6}
        (2, 2),
        (2, 4),
        (2, 6),
        (2, 12),
        # Groups of the form Z/3Z x Z/3NZ for N in {1, 2, 3}
        (3, 3),
        (3, 6),
        (3, 9),
        # The group Z/6Z x Z/6Z, which occurs only over Q(sqrt(-3))
        (6, 6)
    ]

    max_order = 0
    largest_group_factors = None

    print("Analyzing possible non-cyclic torsion subgroups over Q(sqrt(-3)):")
    print("-" * 60)

    for m, n in possible_groups:
        order = m * n
        print(f"Group: Z/{m}Z x Z/{n}Z,  Order Calculation: {m} * {n} = {order}")
        if order > max_order:
            max_order = order
            largest_group_factors = (m, n)
            
    if largest_group_factors is None:
        print("\nNo non-cyclic groups found in the list.")
        # Exit if for some reason the list was empty
        return
        
    print("-" * 60)
    print(f"\nThe largest order of a non-cyclic torsion subgroup is {max_order}.")
    
    m, n = largest_group_factors
    print(f"This maximum order comes from the group Z/{m}Z x Z/{n}Z.")
    print("The final equation for the largest order is:")
    print(f"{m} * {n} = {max_order}")

solve_torsion_order()