import math

def is_prime(n):
    """A simple primality test."""
    if n < 2:
        return False
    for i in range(2, int(math.sqrt(n)) + 1):
        if n % i == 0:
            return False
    return True

def find_next_prime(n):
    """Finds the smallest prime number larger than n."""
    p = n + 1
    while True:
        if is_prime(p):
            return p
        p += 1

def inverse(a, m):
    """Calculates the modular multiplicative inverse of a modulo m."""
    # Using Fermat's Little Theorem: a^(p-2) = a^(-1) (mod p)
    # This is efficient for prime moduli.
    return pow(a, m - 2, m)

def construct_rigid_matrix(N):
    """
    Constructs an N x N Cauchy matrix, a candidate for a rigid matrix.

    Args:
        N: The dimension of the matrix.

    Returns:
        A 2D list representing the N x N Cauchy matrix.
    """
    if N <= 0:
        return []

    # A Cauchy matrix C_ij = 1 / (x_i - y_j) requires 2N distinct elements.
    # We construct it over a finite field F_p where p is a prime > 2*N.
    p = find_next_prime(2 * N)
    print(f"Constructing Cauchy matrix of size {N}x{N} over the finite field F_{p}")

    # Choose two disjoint sets of N distinct elements {x_i} and {y_j}
    x = [i + 1 for i in range(N)]  # x = {1, 2, ..., N}
    y = [i + 1 + N for i in range(N)] # y = {N+1, N+2, ..., 2N}

    matrix = [[0] * N for _ in range(N)]

    print("\nConstructing matrix M where M_ij = (x_i - y_j)^-1 mod p...")
    for i in range(N):
        for j in range(N):
            # Calculate the difference x_i - y_j
            diff = (x[i] - y[j]) % p
            # Calculate the modular inverse
            inv_diff = inverse(diff, p)
            matrix[i][j] = inv_diff
    
    return matrix

def main():
    try:
        N_str = input("Enter the integer dimension N for the matrix: ")
        N = int(N_str)
        if N <= 0:
            raise ValueError("N must be a positive integer.")
    except (ValueError, TypeError) as e:
        print(f"Invalid input: {e}")
        return

    matrix = construct_rigid_matrix(N)

    print("\nConstructed Matrix:")
    for row in matrix:
        print(" ".join(map(str, row)))

if __name__ == "__main__":
    main()
