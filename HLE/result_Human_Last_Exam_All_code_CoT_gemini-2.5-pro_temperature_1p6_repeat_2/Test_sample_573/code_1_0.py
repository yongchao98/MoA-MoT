import math

def calculate_simplices():
    """
    Calculates the number of n-simplices for the given problem parameters.
    
    The number of n-simplices in N_.(Z_N)_{k/} is the number of sequences
    (x_0, ..., x_n) such that -N <= x_0 <= ... <= x_n <= k.
    This is a combinatorial problem of choosing n+1 items from N+k+1 options
    with replacement, which is given by the formula C(N+k+n+1, n+1).
    """
    N = 200
    k = 13

    print(f"Calculating the number of n-simplices for N={N}, k={k}, and n <= 5.")
    print("-" * 60)

    for n in range(6):
        # Parameters for the combination formula C(n, k)
        # We are choosing n+1 numbers from a set of N+k+1 integers.
        # The number of ways is C((N+k+1) + (n+1) - 1, n+1)
        total_items = N + k + n + 1
        items_to_choose = n + 1
        
        # Calculate the number of simplices using math.comb
        result = math.comb(total_items, items_to_choose)
        
        # Print the detailed equation as requested
        print(f"For n = {n}, the number of {n}-simplices is:")
        print(f"C({N} + {k} + {n} + 1, {n} + 1) = C({total_items}, {items_to_choose}) = {result}")
        print()

calculate_simplices()