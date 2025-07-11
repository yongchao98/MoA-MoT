import math

def find_largest_torsion_order():
    """
    This function finds the largest order of a non-cyclic torsion subgroup
    of an elliptic curve over the number field Q(sqrt(-3)).

    The data is based on the complete classification of torsion structures
    for elliptic curves over quadratic fields.
    """

    # List of possible non-cyclic torsion groups over Q(sqrt(-3))
    # Each tuple (n1, n2) represents the group Z/n1Z x Z/n2Z
    # This list is derived from the known classification theorems.
    # Note: n2 must divide n1. The representation here is for convenience.
    
    # We use tuples (a, b) to represent Z/aZ x Z/bZ
    possible_groups = [
        # Groups that can occur over Q
        (2, 2),  # Z/2Z x Z/2Z
        (2, 4),  # Z/2Z x Z/4Z
        (2, 6),  # Z/2Z x Z/6Z
        (2, 8),  # Z/2Z x Z/8Z

        # Other groups that can occur over various quadratic fields, including Q(sqrt(-3))
        (2, 10), # Z/2Z x Z/10Z
        (2, 12), # Z/2Z x Z/12Z
        (3, 3),  # Z/3Z x Z/3Z
        (3, 6),  # Z/3Z x Z/6Z

        # Groups that, among quadratic fields, occur only over Q(sqrt(-3))
        (3, 9),  # Z/3Z x Z/9Z
        (6, 6)   # Z/6Z x Z/6Z
    ]

    max_order = 0
    best_group = None

    print("Analyzing possible non-cyclic torsion groups over Q(sqrt(-3)):\n")

    for group in possible_groups:
        n1, n2 = group
        order = n1 * n2
        
        # Following the prompt's request to output each number in the final equation
        print(f"Group Z/{n1}Z x Z/{n2}Z has order {n1} * {n2} = {order}")

        if order > max_order:
            max_order = order
            best_group = group
    
    print("\n------------------------------------------------------------")
    print("The group with the largest order is Z/{}Z x Z/{}Z.".format(best_group[0], best_group[1]))
    print("The largest possible order is {} * {} = {}.".format(best_group[0], best_group[1], max_order))
    print("------------------------------------------------------------\n")
    
    return max_order

# Execute the function and print the final answer
largest_order = find_largest_torsion_order()
print(f"The largest order of a non-cyclic torsion subgroup of an elliptic curve over Q(sqrt(-3)) is:")
print(f"<<<{largest_order}>>>")
