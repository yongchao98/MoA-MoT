import math

def calculate_rigid_rank(N, delta):
    """
    Calculates the largest rank r for which an FNP algorithm can construct
    a (delta, r)-rigid N x N matrix.

    Args:
        N (int): The dimension of the matrix.
        delta (float): The fraction of entries that can be changed (a small constant).

    Returns:
        None. Prints the result.
    """
    if not (0 < delta < 1):
        print("Error: delta must be a value between 0 and 1.")
        return
    if not (isinstance(N, int) and N > 0):
        print("Error: N must be a positive integer.")
        return

    # The largest rank r is derived from the existence proof: r < N * (1 - sqrt(delta))
    # We take the floor to get the largest integer rank.
    r_float = N * (1 - math.sqrt(delta))
    r = math.floor(r_float)

    # Print the equation and the result
    print("Based on the existence proof for rigid matrices, an FNP algorithm can construct a matrix for a rank 'r' up to:")
    print(f"r = floor(N * (1 - sqrt(delta)))")
    print(f"r = floor({N} * (1 - sqrt({delta})))")
    print(f"r = floor({N} * (1 - {math.sqrt(delta):.4f}))")
    print(f"r = floor({N} * {1 - math.sqrt(delta):.4f})")
    print(f"r = floor({r_float:.4f})")
    print(f"So, the largest integer rank is r = {r}")

# Example usage with N=1000 and a small delta
N = 1000
delta = 0.01
calculate_rigid_rank(N, delta)