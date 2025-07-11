import numpy as np

def solve_and_find_rank():
    """
    This function demonstrates the solution by:
    1. Defining example A, b, and x.
    2. Calculating the perturbation matrix E with the minimum Frobenius norm.
    3. Computing and printing the rank of E.
    """
    # Define matrix A (3x2), vector b (3x1), and vector x (2x1)
    # We choose them such that b - Ax is not zero.
    A = np.array([[1.0, 2.0],
                    [3.0, 4.0],
                    [5.0, 6.0]])

    b = np.array([[1.0],
                    [1.0],
                    [1.0]])

    x = np.array([[1.0],
                    [1.0]])

    print("Given matrices:")
    print(f"A =\n{A}")
    print(f"b =\n{b}")
    print(f"x =\n{x}")
    print("-" * 30)

    # Calculate r = b - Ax
    r = b - A @ x

    # Calculate the squared L2-norm of x
    x_norm_sq = (np.linalg.norm(x)**2)

    # Calculate the transpose of x
    x_transpose = x.T

    # Calculate E = r * x^T / ||x||^2
    E = (r @ x_transpose) / x_norm_sq
    
    # Calculate the rank of E
    rank_E = np.linalg.matrix_rank(E)

    print("Calculation of E based on E = r * x_transpose / ||x||^2:")
    print(f"r = b - Ax =\n{r}")
    print(f"x_transpose =\n{x_transpose}")
    print(f"x_norm_sq = ||x||^2 = {x_norm_sq:.4f}")
    print("-" * 30)
    
    print(f"The resulting matrix E is:\n{E}")
    print("-" * 30)
    
    print(f"The rank of E is: {rank_E}")

solve_and_find_rank()

# The greatest possible rank of E is 1.
# We have shown a case where the rank is 1.
# The rank can also be 0 if r = b - Ax is the zero vector.
# Thus, the maximum possible rank is 1.