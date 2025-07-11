def calculate_rank():
    """
    Calculates the rank of the kernel of the Hurewicz homomorphism for the given space Y.
    The formula for the rank r of the commutator subgroup of a free product of n finite groups G_i is:
    r = 1 - sum(product(|G_j|) for j != i) + (n-1) * product(|G_i|)
    """
    # Orders of the fundamental groups of X1, X2, X3
    g1_order = 5
    g2_order = 8
    g3_order = 2
    n = 3

    # Product of all orders
    prod_all = g1_order * g2_order * g3_order

    # Sum of products of orders, leaving one group out each time
    # Equivalent to (prod_all / g1_order) + (prod_all / g2_order) + (prod_all / g3_order)
    sum_of_prods = (g2_order * g3_order) + (g1_order * g3_order) + (g1_order * g2_order)

    # The rank formula
    rank = 1 - sum_of_prods + (n - 1) * prod_all

    # Print the calculation steps
    print("The rank is calculated using the formula:")
    print("Rank = 1 - (sum of products of orders of group combinations) + (n-1) * (product of all group orders)")
    print(f"Orders of the groups are |G1|={g1_order}, |G2|={g2_order}, |G3|={g3_order}.")
    print(f"The number of groups is n={n}.")
    print("\nStep 1: Calculate the product of all group orders:")
    print(f"Product = {g1_order} * {g2_order} * {g3_order} = {prod_all}")
    
    print("\nStep 2: Calculate the sum of products of orders, leaving one group out:")
    calc_sum_str = f"({g2_order} * {g3_order}) + ({g1_order} * {g3_order}) + ({g1_order} * {g2_order}) = {g2_order*g3_order} + {g1_order*g3_order} + {g1_order*g2_order}"
    print(f"Sum = {calc_sum_str} = {sum_of_prods}")

    print("\nStep 3: Apply the full formula:")
    print(f"Rank = 1 - {sum_of_prods} + ({n} - 1) * {prod_all}")
    print(f"Rank = 1 - {sum_of_prods} + {n-1} * {prod_all}")
    print(f"Rank = 1 - {sum_of_prods} + {(n - 1) * prod_all}")
    print(f"Rank = {1 - sum_of_prods} + {(n - 1) * prod_all}")
    print(f"Rank = {rank}")

calculate_rank()