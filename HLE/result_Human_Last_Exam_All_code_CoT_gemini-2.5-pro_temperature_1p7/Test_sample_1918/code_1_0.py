import math

def max_complementary_choi_rank(d, n, r):
    """
    Calculates the maximal rank of the Choi matrix of the complementary channel.

    Args:
        d (int): Dimension of the input Hilbert space H1.
        n (int): Dimension of the output Hilbert space H2.
        r (int): Rank of the Choi matrix of the channel Lambda.
    
    Returns:
        int: The maximal rank of the Choi matrix of the complementary channel Lambda^c.
    """
    # A quantum channel with these parameters is only possible if |d-r| <= n.
    # We assume the provided (d, n, r) triplet is valid.
    if d <= 0 or n <= 0 or r <= 0:
        raise ValueError("Dimensions and rank must be positive integers.")
        
    # The maximal rank is given by the formula min(n, d + r - 1)
    result = min(n, d + r - 1)
    
    print(f"Given a quantum channel from a Hilbert space of dimension d={d} to one of dimension n={n},")
    print(f"with its Choi matrix having a rank of r={r},")
    print("the maximal rank of the Choi matrix of the complementary channel is calculated as:")
    print(f"max_rank = min(n, d + r - 1)")
    print(f"max_rank = min({n}, {d} + {r} - 1)")
    print(f"max_rank = min({n}, {d + r - 1})")
    print(f"max_rank = {result}")

# Example usage with some arbitrary valid parameters.
# Case 1: n is the smaller value
print("--- Case 1 ---")
max_complementary_choi_rank(d=10, n=8, r=5)

print("\n--- Case 2 ---")
# Case 2: d + r - 1 is the smaller value
max_complementary_choi_rank(d=4, n=20, r=3)
