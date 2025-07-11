import math

def calculate_rigid_rank(N, delta):
    """
    Calculates the largest rank r for which an FNP algorithm can construct
    a (delta, r)-rigid N x N matrix over F_2, based on a known existence proof.

    The existence of such matrices is proven by a counting argument. The result
    of this argument is that for a matrix to be (delta, r)-rigid, the rank r
    can be as large as c*N, where c is a constant dependent on delta.
    This function calculates r for given N and delta.

    Args:
        N (int): The dimension of the square matrix.
        delta (float): The fraction of entries that can be changed.
                       Must be in the range (0, 0.5).

    Returns:
        int: The calculated rank r, or an error message string.
    """
    if not (0 < delta < 0.5):
        return "Error: delta must be between 0 and 0.5 for this formula to yield a positive rank."

    # Step 1: Calculate the binary entropy H(delta).
    # H(delta) = -delta*log2(delta) - (1-delta)*log2(1-delta)
    h_delta = -delta * math.log2(delta) - (1 - delta) * math.log2(1 - delta)

    # Step 2: Calculate the constant c for the rank r = c*N.
    # The counting argument guarantees existence for c < 1 - sqrt(H(delta)).
    c = 1 - math.sqrt(h_delta)
    
    # Step 3: Calculate the rank r.
    r = math.floor(c * N)

    print(f"For N = {N} and delta = {delta}:")
    print("The existence of rigid matrices is guaranteed by a counting argument.")
    print("The condition for existence over F_2 is r/N < 1 - sqrt(H(delta)), where H(delta) is the binary entropy.")
    print(f"1. Binary Entropy H({delta}) = {h_delta:.4f}")
    print(f"2. Constant c = 1 - sqrt(H({delta})) = 1 - sqrt({h_delta:.4f}) = {c:.4f}")
    print(f"3. Final Equation for rank r: r = floor(c * N) = floor({c:.4f} * {N})")
    print(f"The resulting rank is: r = {r}")
    
    return r

if __name__ == '__main__':
    # Example usage of the function.
    # The user can modify these values.
    N_val = 1000
    delta_val = 0.01

    calculate_rigid_rank(N_val, delta_val)