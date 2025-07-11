import math

def solve_polynomial_complexity():
    """
    Calculates the smallest complexity of a matrix product for the polynomial f_n,d.

    The user can modify the values of n and d below to get the complexity
    for a different case.
    """
    # Per the problem description, n and d are integers with 2 <= d <= n and d is odd.
    # Example values:
    # n=4, d=3
    # n=5, d=3
    # n=5, d=5
    # To run this code for a different case, please modify the n and d variables below.
    n = 4
    d = 3
    
    if not (isinstance(n, int) and isinstance(d, int) and 2 <= d <= n and d % 2 != 0):
        print(f"Error: Input values must satisfy the constraints: n and d are integers, 2 <= d <= n, and d is odd.")
        print(f"Current values are: n={n}, d={d}")
        return

    def combinations(n, k):
        """
        Calculates the binomial coefficient C(n, k) using integer arithmetic.
        """
        # This function is equivalent to math.comb(n, k) available in Python 3.8+
        if k < 0 or k > n:
            return 0
        if k == 0 or k == n:
            return 1
        if k > n // 2:
            k = n - k
        
        res = 1
        for i in range(k):
            res = res * (n - i) // (i + 1)
        return res

    # The formula for complexity is based on a symmetric matrix construction.
    # Let k = (d-1)/2. The complexity is 2 + 2 * sum_{i=1 to k} C(n,i).
    k = (d - 1) // 2

    if k == 0:
        # This case happens if d=1, but the problem states d >= 2.
        # If it were to happen, the sum would be empty (equals 0).
        complexity = 2
        print(f"2 + 2 * 0 = 2")
        return
    
    # Calculate the sum of binomial coefficients
    coefficients = [combinations(n, i) for i in range(1, k + 1)]
    sum_of_coeffs = sum(coefficients)
    
    # Calculate the final complexity
    complexity = 2 + 2 * sum_of_coeffs
    
    # Format the output string to show all the numbers in the final equation.
    coeffs_str_list = [str(c) for c in coefficients]
    sum_str = " + ".join(coeffs_str_list)

    if k == 1:
        # For a single term, omit the parentheses around the sum for clarity
        print(f"2 + 2 * {sum_str} = 2 + 2 * {sum_of_coeffs} = {complexity}")
    else:
        print(f"2 + 2 * ({sum_str}) = 2 + 2 * {sum_of_coeffs} = {complexity}")

solve_polynomial_complexity()
<<<10>>>