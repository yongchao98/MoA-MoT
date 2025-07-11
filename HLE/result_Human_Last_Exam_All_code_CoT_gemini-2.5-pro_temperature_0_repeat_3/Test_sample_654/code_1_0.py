import math

def solve_polynomial_coeffs():
    """
    This function solves the problem for given p, k, n.
    Since the problem asks for the answer in terms of p, k, and n,
    we will ask the user for these inputs.
    """
    try:
        p_str = input("Enter the odd prime p: ")
        p = int(p_str)
        if p < 3 or any(p % i == 0 for i in range(2, int(math.sqrt(p)) + 1)):
            print("p must be an odd prime.")
            return

        k_str = input("Enter the integer k (k >= 1): ")
        k = int(k_str)
        if k < 1:
            print("k must be >= 1.")
            return

        n_str = input("Enter the integer n (n >= 1): ")
        n = int(n_str)
        if n < 1:
            print("n must be >= 1.")
            return

    except ValueError:
        print("Invalid input. Please enter integers.")
        return

    # The number of coefficients not divisible by p^k is given by the formula:
    # p^n + 1 if n < k
    # p^(k-1) + 1 if n >= k
    # This can be written compactly as p^(min(n, k-1)) + 1

    if n < k:
        exponent = n
        result = p**exponent + 1
        print(f"Since n < k ({n} < {k}), the formula is p^n + 1.")
        print(f"The number of coefficients is {p}^{n} + 1 = {result}")
    else: # n >= k
        exponent = k - 1
        result = p**exponent + 1
        print(f"Since n >= k ({n} >= {k}), the formula is p^(k-1) + 1.")
        print(f"The number of coefficients is {p}^{k-1} + 1 = {result}")

solve_polynomial_coeffs()