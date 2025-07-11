def calculate_kernel_rank():
    """
    Calculates the rank of the kernel of the Hurewicz homomorphism for the given space Y.
    Y is the connected sum of X1, X2, and X3.
    pi_1(X1) = Z_5, pi_1(X2) = Z_8, pi_1(X3) = Z_2.
    The kernel K is the kernel of the map G -> G^ab, where G = Z_5 * Z_8 * Z_2.
    The rank is computed using a formula for the kernel of the map from a free product
    to a direct product of finite groups.
    """

    # Orders of the cyclic groups
    n1 = 5
    n2 = 8
    n3 = 2
    groups_orders = [n1, n2, n3]

    # Number of groups in the free product
    n = len(groups_orders)

    # Order of the direct product H = Z_5 x Z_8 x Z_2
    H_order = 1
    for order in groups_orders:
        H_order *= order

    # Sum of the reciprocals of the group orders
    sum_inv_orders = sum(1/order for order in groups_orders)

    # The formula for the rank of the kernel K
    # rank(K) = 1 + |H| * ( (n-1) - sum(1/|G_i|) )
    rank = 1 + H_order * ((n - 1) - sum_inv_orders)

    # Output the final equation with numbers filled in
    print(f"The rank is calculated using the formula: rank = 1 + |H| * ( (n-1) - sum(1/|G_i|) )")
    print(f"Plugging in the numbers:")
    # The string construction shows the full equation as requested
    equation_str = f"rank = 1 + {H_order} * (({n} - 1) - (1/{n1} + 1/{n2} + 1/{n3}))"
    print(equation_str)

    # Perform the intermediate calculation steps and print them
    step1 = f"rank = 1 + {H_order} * ({n-1} - {sum_inv_orders})"
    print(step1)
    
    step2_val = (n-1) - sum_inv_orders
    step2 = f"rank = 1 + {H_order} * ({step2_val})"
    print(step2)

    step3_val = H_order * step2_val
    step3 = f"rank = 1 + {step3_val}"
    print(step3)
    
    # Final result
    print("\nThe final rank is:")
    print(int(rank))

# Run the calculation
calculate_kernel_rank()