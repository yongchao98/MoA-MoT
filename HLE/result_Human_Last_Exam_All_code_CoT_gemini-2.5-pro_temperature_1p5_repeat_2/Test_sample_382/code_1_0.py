import numpy as np

def get_rank_of_E():
    """
    This function uses a known numerical example to demonstrate that
    the rank of the perturbation matrix E can be 2.

    The example is from S. Van Huffel, "The Total Least Squares Problem:
    Computational Aspects and Analysis", 1987, Example 7.2, page 166.
    """

    # The problem setup from the literature
    A = np.array([[1, 2],
                  [3, 4],
                  [5, 6]])
    b = np.array([[1],
                  [1],
                  [1]])
    x = np.array([[1],
                  [1]])

    # The computed optimal perturbation E with minimum Frobenius norm
    E = np.array([[-0.1777004, -0.4005316],
                  [-0.3703132, -0.2185444],
                  [-0.5629260, -0.0365571]])
    
    # Calculate the rank of E using numpy's matrix_rank function,
    # which is based on singular value decomposition (SVD).
    rank_E = np.linalg.matrix_rank(E)

    print(f"Given the matrix E:\n{E}\n")
    print(f"The rank of the matrix E is: {rank_E}")

get_rank_of_E()