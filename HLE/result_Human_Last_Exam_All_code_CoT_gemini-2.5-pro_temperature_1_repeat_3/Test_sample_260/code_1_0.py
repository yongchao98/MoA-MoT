import math

def solve_rank():
    """
    Calculates the rank of the kernel K based on the derived formula.
    """
    # Orders of the cyclic groups corresponding to X1, X2, X3
    orders = [5, 8, 2]
    k = len(orders)

    # Calculate the order of the abelianization H_1(Y)
    A_abs = 1
    for n in orders:
        A_abs *= n

    # Calculate the sum of the reciprocals of the orders
    sum_inv_orders = sum(1/n for n in orders)

    # Apply the formula for the rank
    # rank = 1 + |A| * ( (k-1) - sum(1/|G_i|) )
    rank = 1 + A_abs * (k - 1 - sum_inv_orders)

    # We can perform the calculation with integers to avoid floating point issues.
    # The sum of inverses is (1/5 + 1/8 + 1/2) = (8 + 5 + 20)/40 = 33/40.
    # rank = 1 + 80 * (2 - 33/40) = 1 + 80 * ((80-33)/40) = 1 + 80 * (47/40)
    # rank = 1 + 2 * 47 = 1 + 94 = 95
    
    print(f"The number of spaces is k = {k}")
    print(f"The orders of the fundamental groups are {orders}")
    print(f"The order of H_1(Y) is |A| = 5 * 8 * 2 = {A_abs}")
    print(f"The sum of the reciprocal orders is 1/5 + 1/8 + 1/2 = {sum_inv_orders}")
    print(f"The rank is calculated using the formula: rank = 1 + |A| * (k - 1 - sum_inv_orders)")
    print(f"rank = 1 + {A_abs} * (({k} - 1) - (1/{orders[0]} + 1/{orders[1]} + 1/{orders[2]}))")
    print(f"rank = 1 + {A_abs} * ({k-1} - {sum_inv_orders})")
    print(f"rank = 1 + {A_abs} * {k - 1 - sum_inv_orders}")
    print(f"rank = {rank}")

    # The result should be an integer
    final_rank = int(round(rank))
    print(f"\nThe rank of the kernel K as a free group is: {final_rank}")

solve_rank()