import sys

def matrix_multiply(A, B, m):
    """Multiplies two 2x2 matrices modulo m."""
    C = [[0, 0], [0, 0]]
    C[0][0] = (A[0][0] * B[0][0] + A[0][1] * B[1][0]) % m
    C[0][1] = (A[0][0] * B[0][1] + A[0][1] * B[1][1]) % m
    C[1][0] = (A[1][0] * B[0][0] + A[1][1] * B[1][0]) % m
    C[1][1] = (A[1][0] * B[0][1] + A[1][1] * B[1][1]) % m
    return C

def matrix_power(A, n, m):
    """Computes A^n modulo m for a 2x2 matrix A using exponentiation by squaring."""
    res = [[1, 0], [0, 1]]  # Identity matrix
    base = A
    while n > 0:
        if n % 2 == 1:
            res = matrix_multiply(res, base, m)
        base = matrix_multiply(base, base, m)
        n //= 2
    return res

def fib_mod(n, m):
    """Computes the n-th Fibonacci number modulo m."""
    if m == 0:
        # This case is not well-defined for this problem, but handle it.
        raise ValueError("Modulo cannot be zero.")
    if n < 0:
        raise ValueError("Fibonacci index must be non-negative.")
    if n == 0:
        return 0
    
    # The Fibonacci sequence can be generated using matrix exponentiation:
    # M = [[1, 1], [1, 0]]
    # M^k = [[F(k+1), F(k)], [F(k), F(k-1)]]
    # So, F(n) is the top-right element of M^(n).
    M = [[1, 1], [1, 0]]
    M_pow_n = matrix_power(M, n, m)
    return M_pow_n[0][1]

def find_next_term(start_k):
    """
    Finds the next integer k > start_k such that k divides L(k+1).
    L(n) = 2*F(n+1) - 1, so we check for (2*F(k+2) - 1) % k == 0.
    """
    k = start_k
    while True:
        k += 1
        # To avoid slow execution, we can add a progress indicator,
        # but the search is expected to be relatively fast.
        # Let's check about 15000 numbers which should be fast enough
        if k > start_k + 15000:
             print("Search taking longer than expected, no solution found in the near range.")
             return None

        # We need F(k+2) mod k
        f_k_plus_2 = fib_mod(k + 2, k)
        
        # Condition is L(k+1) % k == 0, which is (2 * F(k+2) - 1) % k == 0
        if (2 * f_k_plus_2 - 1) % k == 0:
            return k

# The original sequence
sequence = [24663, 35005, 119261, 196219, 211770, 227296]
last_term = sequence[-1]

print(f"Searching for the next integer k > {last_term} that satisfies the condition...")

next_term = find_next_term(last_term)

if next_term:
    sequence.append(next_term)
    print("\nFound the next term. The completed sequence is:")
    for number in sequence:
        print(number)
