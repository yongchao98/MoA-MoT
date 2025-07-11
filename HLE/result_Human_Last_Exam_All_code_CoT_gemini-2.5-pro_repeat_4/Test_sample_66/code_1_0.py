import numpy as np

def create_tonga_matrix(M, N, m, n, a, b):
    """
    Creates a matrix representing the Tongan flag with maximal rank geometry.

    Args:
        M (int): Total height of the flag matrix.
        N (int): Total width of the flag matrix.
        m (int): Height of the canton.
        n (int): Width of the canton.
        a (float): Value for red pixels.
        b (float): Value for white pixels.

    Returns:
        np.ndarray: The matrix representing the flag.
        np.ndarray: The W matrix (white parts of the canton).
    """
    # Start with a full red flag
    flag_matrix = np.full((M, N), a)

    # Create the white canton
    flag_matrix[:m, :n] = b

    # Create the red cross in the canton.
    # To maximize rank, we use a cross shape made of two partially
    # overlapping rectangles H and V.
    
    # Define vertical bar V region
    i_v1, i_v2 = m // 5, 3 * m // 5
    j_v1, j_v2 = n // 5, 3 * n // 5
    
    # Define horizontal bar H region
    i_h1, i_h2 = 2 * m // 5, 4 * m // 5
    j_h1, j_h2 = 2 * n // 5, 4 * n // 5

    # Create a boolean mask for the cross
    cross_mask = np.zeros((m, n), dtype=bool)
    # Mark horizontal bar region
    cross_mask[i_h1:i_h2, j_h1:j_h2] = True
    # Mark vertical bar region
    cross_mask[i_v1:i_v2, j_v1:j_v2] = True

    # Apply the cross to the canton
    flag_matrix[:m, :n][cross_mask] = a

    # Create the W matrix (1 for white, 0 for red in canton)
    canton = flag_matrix[:m, :n]
    W = (canton == b).astype(int)
    
    return flag_matrix, W

def solve():
    """
    Calculates and prints the maximal rank of the Tongan flag matrix.
    """
    # Define dimensions (must be large enough for the geometry)
    M, N = 100, 150  # Flag dimensions
    m, n = 50, 75    # Canton dimensions

    # Choose a and b to maximize rank (a!=b, a!=0)
    a, b = 1.0, 2.0

    # Create the matrices
    A, W = create_tonga_matrix(M, N, m, n, a, b)

    # Calculate the ranks
    rank_A = np.linalg.matrix_rank(A)
    rank_W = np.linalg.matrix_rank(W)

    # Print the explanation and the result
    print("The rank of the flag matrix 'A' can be decomposed as rank(A) = rank(W) + 1,")
    print("where W is the matrix representing the white areas of the canton (1 for white, 0 for red).")
    print("To maximize the rank of A, we need to maximize the rank of W.")
    print("By choosing a cross geometry that creates four distinct row/column patterns, we can achieve rank(W) = 4.")
    print("\nCalculating the ranks with this optimized geometry:")
    
    # Output the numbers in the final equation
    print(f"Rank of the canton's white-pattern matrix, rank(W) = {rank_W}")
    print(f"Maximal rank of the flag matrix, rank(A) = rank(W) + 1 = {rank_W} + 1 = {rank_A}")

solve()
<<<5>>>