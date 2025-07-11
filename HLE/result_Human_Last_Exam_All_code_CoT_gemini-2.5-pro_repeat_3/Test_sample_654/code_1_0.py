import math

def solve_polynomial_coefficients(p, k, n):
    """
    Computes the number of coefficients not divisible by p^k in the given polynomial sequence.

    Args:
        p (int): An odd prime number (p >= 3).
        k (int): An integer (k >= 1).
        n (int): An integer (n >= 1).

    Returns:
        int: The number of coefficients not divisible by p^k.
    """
    # Verify constraints from the problem statement
    if not (p >= 3 and p % 2 == 1):
        # A simple primality test for demonstration. For larger p, a more robust test would be needed.
        is_prime = all(p % i != 0 for i in range(3, int(math.sqrt(p)) + 1, 2))
        if not is_prime and p != 3:
             print(f"Warning: p={p} may not be an odd prime.")

    if not (k >= 1):
        raise ValueError("k must be an integer >= 1")
    if not (n >= 1):
        raise ValueError("n must be an integer >= 1")

    # Based on the derived case analysis:
    if k == 1:
        # Case 1: k=1
        # The number of coefficients is 2.
        result = 2
        print(f"For p={p}, k={k}, n={n}, the case is k=1.")
        print(f"The number of coefficients not divisible by p^k is 2.")
        
    elif n < k:
        # Case 2: k > 1 and n < k
        # The number of coefficients is p^n + 1.
        result = p**n + 1
        print(f"For p={p}, k={k}, n={n}, the case is k > 1 and n < k.")
        print(f"The result is p^n + 1 = {p}^{n} + 1 = {result}.")

    else: # k > 1 and n >= k
        # Case 3: k > 1 and n >= k
        # The number of coefficients is p^(k-1) + 1.
        result = p**(k - 1) + 1
        print(f"For p={p}, k={k}, n={n}, the case is k > 1 and n >= k.")
        print(f"The result is p^(k-1) + 1 = {p}^({k-1}) + 1 = {result}.")
        
    return result

# You can change these values to test with different parameters.
p_val = 5
k_val = 3
n_val = 4

# Execute the solution with example values.
solve_polynomial_coefficients(p_val, k_val, n_val)