def solve_rank_of_kernel():
    """
    Calculates the rank of the kernel of the Hurewicz homomorphism for the given space Y.

    The space Y is the connected sum of three spaces X1, X2, X3 with fundamental groups
    Z_5, Z_8, and Z_2, respectively.
    The fundamental group of Y is the free product G = Z_5 * Z_8 * Z_2.
    The kernel K of the Hurewicz map is the commutator subgroup [G, G], which is a free group.
    Its rank r can be calculated with the formula:
    r = 1 + (k-1) * Product(n_i) - Sum(Product(n_i) / n_j)
    where n_i are the orders of the cyclic groups and k is the number of groups.
    """
    
    # Orders of the cyclic groups corresponding to pi_1(X_1), pi_1(X_2), and pi_1(X_3)
    orders = [5, 8, 2]
    n1, n2, n3 = orders
    
    # Number of groups being combined in the free product
    k = len(orders)
    
    # Product of the orders
    prod_n = 1
    for n in orders:
        prod_n *= n
        
    # Summation term in the formula
    sum_term = 0
    for n in orders:
        sum_term += prod_n // n
        
    # Calculate the rank using the formula
    rank = 1 + (k - 1) * prod_n - sum_term
    
    # Print out the detailed calculation
    print("The fundamental groups of X1, X2, and X3 are Z_5, Z_8, and Z_2.")
    print(f"The group G = pi_1(Y) is the free product of these {k} groups.")
    print("The kernel K is the commutator subgroup [G, G], which is a free group.")
    print("The rank 'r' is calculated using the formula: r = 1 + (k-1)*P - S")
    print(f"where P is the product of the orders {orders}, and S is the sum of P divided by each order.")
    
    print("\n--- Calculation ---")
    print(f"Orders (n1, n2, n3): ({n1}, {n2}, {n3})")
    print(f"Number of groups (k): {k}")
    print(f"Product of orders (P = n1*n2*n3): {n1} * {n2} * {n3} = {prod_n}")
    
    term1_val = prod_n // n1
    term2_val = prod_n // n2
    term3_val = prod_n // n3
    print(f"Sum term (S = P/n1 + P/n2 + P/n3): {term1_val} + {term2_val} + {term3_val} = {sum_term}")
    
    k_minus_1_prod = (k - 1) * prod_n
    print(f"Rank = 1 + (k-1)*P - S = 1 + ({k}-1)*{prod_n} - {sum_term}")
    print(f"     = 1 + {k_minus_1_prod} - {sum_term}")
    print(f"     = {1 + k_minus_1_prod} - {sum_term}")
    print(f"     = {rank}")
    print("\nThe rank of the kernel K as a free group is:")
    print(rank)

solve_rank_of_kernel()