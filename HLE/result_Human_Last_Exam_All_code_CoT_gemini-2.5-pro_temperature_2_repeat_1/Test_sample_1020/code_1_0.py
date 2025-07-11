import math

def get_smallest_complexity(n, d):
    """
    Calculates the smallest complexity of a matrix product computing f_n,d.

    Args:
        n: The upper bound for the variable indices.
        d: The number of variables in each term of the polynomial.
    """
    
    # Input validation as per the problem description
    if not (isinstance(n, int) and isinstance(d, int)):
        print("Error: n and d must be integers.")
        return
    if not (2 <= d <= n):
        print(f"Error: Constraints not met. We must have 2 <= d <= n, but got d={d} and n={n}.")
        return
    if d % 2 == 0:
        print(f"Error: d must be odd, but got d={d}.")
        return

    def combinations(n, k):
        """Helper function to calculate binomial coefficients C(n,k)"""
        if k < 0 or k > n:
            return 0
        if k == 0 or k == n:
            return 1
        # Use the fact that C(n,k) = C(n, n-k) to optimize
        if k > n // 2:
            k = n - k
        
        # Calculate C(n,k) iteratively to avoid large intermediate numbers
        res = 1
        for i in range(k):
            res = res * (n - i) // (i + 1)
        return res

    # As d is odd, let d = 2k+1. Then k = (d-1)/2.
    k = (d - 1) // 2
    
    sum_of_combs = 0
    # The list to hold the string representation of each C(n,i) value
    comb_values_str_list = []

    for i in range(1, k + 1):
        comb_val = combinations(n, i)
        sum_of_combs += comb_val
        comb_values_str_list.append(str(comb_val))

    total_complexity = 2 + 2 * sum_of_combs
    
    # Building the equation string with evaluated numbers
    equation_part = " + ".join(comb_values_str_list)
    
    print(f"For n={n} and d={d}, the smallest complexity is calculated as follows:")
    print(f"The formula is: 2 + 2 * (C(n,1) + ... + C(n,k)) where k = (d-1)/2.")
    print(f"Here, k = ({d}-1)/2 = {k}.")
    print("\nFinal calculation:")
    # Print the equation with all the numbers
    if k > 0:
        print(f"Result = 2 + 2 * ({equation_part})")
        print(f"       = 2 + 2 * {sum_of_combs}")
        print(f"       = 2 + {2 * sum_of_combs}")
        print(f"       = {total_complexity}")
    else: # This happens if d=1, but the problem states d>=2. As a safeguard:
        print(f"Result = {total_complexity}")


# Example values provided by the prompt. Let's use some illustrative ones.
n_val = 10
d_val = 5

get_smallest_complexity(n_val, d_val)

<<<112>>>