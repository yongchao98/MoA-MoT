def solve_rank_of_kernel():
    """
    Calculates the rank of the kernel of the Hurewicz homomorphism for the
    connected sum of the given topological spaces.
    """
    # Orders of the fundamental groups of X1, X2, X3
    n1 = 5
    n2 = 8
    n3 = 2
    orders = [n1, n2, n3]
    
    # Number of spaces being connected-summed
    k = len(orders)
    
    # Calculate the product of the orders
    product_of_orders = 1
    for n in orders:
        product_of_orders *= n
        
    # Calculate the sum of the reciprocals of the orders
    sum_of_reciprocals = 0
    for n in orders:
        sum_of_reciprocals += 1/n
        
    # Apply the formula for the rank of the commutator subgroup
    # rank = 1 + (product of orders) * ( (k-1) - sum of reciprocals )
    rank = 1 + product_of_orders * ( (k - 1) - sum_of_reciprocals )
    
    # The rank must be an integer. We round to handle potential floating point inaccuracies.
    final_rank = int(round(rank))

    # Print the calculation steps
    print("The fundamental groups of the spaces are Z_5, Z_8, and Z_2.")
    print("The fundamental group of their connected sum Y is G = Z_5 * Z_8 * Z_2.")
    print("The kernel K of the Hurewicz map is the commutator subgroup [G, G].")
    print("The rank of K is calculated using the formula:")
    print("rank = 1 + (n1 * n2 * n3) * ((k - 1) - (1/n1 + 1/n2 + 1/n3))")
    print(f"rank = 1 + ({n1} * {n2} * {n3}) * (({k} - 1) - (1/{n1} + 1/{n2} + 1/{n3}))")
    print(f"rank = 1 + {product_of_orders} * ({k-1} - {sum_of_reciprocals})")
    print(f"rank = 1 + {product_of_orders} * ({k-1 - sum_of_reciprocals})")
    print(f"rank = 1 + {product_of_orders * (k-1 - sum_of_reciprocals)}")
    print(f"rank = {final_rank}")
    
    # The final answer is returned separately for the platform
    # print(f"<<<{final_rank}>>>")

solve_rank_of_kernel()