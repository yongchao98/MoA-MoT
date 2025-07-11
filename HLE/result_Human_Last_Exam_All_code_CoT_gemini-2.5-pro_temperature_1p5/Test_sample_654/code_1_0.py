import math

def solve():
    """
    Solves the problem of finding the number of coefficients not divisible by p^k.
    The user is prompted to enter the values of p, k, and n.
    """
    try:
        p_str = input("Enter the odd prime p: ")
        p = int(p_str)
        if p < 3 or p % 2 == 0:
            # A simple primality test for small numbers would be more robust,
            # but we follow the problem's premise that p is an odd prime.
            print("p must be an odd prime (p >= 3).")
            return

        k_str = input("Enter the integer k (k >= 1): ")
        k = int(k_str)
        if k < 1:
            print("k must be a positive integer (k >= 1).")
            return

        n_str = input("Enter the integer n (n >= 1): ")
        n = int(n_str)
        if n < 1:
            # The problem states m>=1, but it is defined as p^n. This implies n>=0.
            # However, the structure of the problem suggests n >= 1. We allow n>=0 for generality.
             print("n must be a non-negative integer (n >= 0). Let's assume n>=1 from context.")
             if n < 0: return

    except ValueError:
        print("Invalid input. Please enter integers only.")
        return

    # The number of coefficients not divisible by p^k is given by the formula:
    # p^min(n, k-1) + 1
    
    # Calculate the exponent
    exponent = min(n, k - 1)
    
    # Calculate the final result
    # We use Python's built-in pow() for integer exponentiation
    result = pow(p, exponent) + 1
    
    # Print the final answer
    print("\n--- Calculation ---")
    print(f"p = {p}")
    print(f"k = {k}")
    print(f"n = {n}")
    print(f"The number of coefficients is given by the formula: p^min(n, k-1) + 1")
    print(f"min(n, k-1) = min({n}, {k}-1) = min({n}, {k-1}) = {exponent}")
    print(f"Result = {p}^{exponent} + 1 = {pow(p, exponent)} + 1 = {result}")
    print("\nFinal Answer:")
    print(result)

solve()