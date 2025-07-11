import math

def calculate_max_rank_complementary_channel(d, n, r):
    """
    Calculates the maximal rank of the Choi matrix of the complementary channel.

    Args:
        d (int): Dimension of the input Hilbert space.
        n (int): Dimension of the output Hilbert space.
        r (int): Rank of the Choi matrix of the original channel.

    Returns:
        int: The maximal rank of the Choi matrix of the complementary channel.
    """
    # A quantum channel from a d-dim to an n-dim space with Choi rank r
    # can only exist if d <= n * r. We check this condition first.
    if d > n * r:
        print(f"Error: A channel with parameters d={d}, n={n}, r={r} cannot exist because d > n*r.")
        return None

    # The maximal rank of the Choi matrix of the complementary channel
    # is given by the formula min(n, d + r).
    max_rank = min(n, d + r)
    return max_rank

# Example values
d_example = 4 # dimension of the input Hilbert space H1
n_example = 5 # dimension of the output Hilbert space H2
r_example = 3 # rank of the Choi matrix of the channel Lambda

# Calculate the result
max_rank_c = calculate_max_rank_complementary_channel(d_example, n_example, r_example)

if max_rank_c is not None:
    # Print the result in a descriptive way showing the equation
    print("The maximal rank of the Choi matrix of the complementary channel is found using the formula:")
    print("max_rank = min(n, d + r)")
    print("\nFor the given example values:")
    print(f"d = {d_example}")
    print(f"n = {n_example}")
    print(f"r = {r_example}")
    print("\nThe calculation is:")
    # We output each number in the final equation as requested
    print(f"max_rank = min({n_example}, {d_example} + {r_example})")
    print(f"max_rank = min({n_example}, {d_example + r_example})")
    print(f"max_rank = {max_rank_c}")
