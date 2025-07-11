import math

def calculate_max_r(N, delta):
    """
    Calculates the largest rank r for which an FNP algorithm can construct
    an (N x N) (delta, r)-rigid matrix.

    Args:
        N (int): The dimension of the matrix.
        delta (float): The fraction of entries that must be changed (0 < delta < 1).
    """
    if not isinstance(N, int) or N <= 0:
        print("Error: N must be a positive integer.")
        return
    if not isinstance(delta, float) or not (0 < delta < 1):
        print("Error: delta must be a float between 0 and 1.")
        return

    print("The formula for the largest rank r is: r = N * (1 - sqrt(delta))")
    print(f"For the given inputs N = {N} and delta = {delta}:")

    # Step 1: Calculate sqrt(delta)
    sqrt_delta = math.sqrt(delta)
    print(f"r = {N} * (1 - sqrt({delta}))")
    print(f"r = {N} * (1 - {sqrt_delta:.4f})")

    # Step 2: Calculate (1 - sqrt(delta))
    one_minus_sqrt_delta = 1 - sqrt_delta
    print(f"r = {N} * {one_minus_sqrt_delta:.4f}")

    # Step 3: Calculate the final value of r
    r = N * one_minus_sqrt_delta
    # Since rank must be an integer, we take the floor.
    final_r = math.floor(r)
    print(f"r = {r:.4f}")
    print(f"\nThe largest integer rank r is: {final_r}")

if __name__ == '__main__':
    try:
        # Example usage with some typical values
        input_N = 1000
        input_delta = 0.01
        calculate_max_r(input_N, input_delta)
    except Exception as e:
        print(f"An error occurred: {e}")
