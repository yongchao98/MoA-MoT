import numpy as np
from scipy.linalg import null_space

def solve_greatest_rank_problem():
    """
    This function demonstrates that the greatest possible rank of the matrix E is 1.
    It constructs a random problem instance (A, b, x), calculates the
    minimum-norm perturbation E such that x is a least-squares solution for (A+E, b),
    and then prints the rank of E.
    """
    # Define problem dimensions
    m, n = 5, 3

    # Create a random problem instance (A, b, x)
    # We ensure that x is not an LS solution for (A,b) to get a non-zero E.
    while True:
        A = np.random.randn(m, n)
        b = np.random.randn(m, 1)
        x = np.random.randn(n, 1)
        
        # Calculate the residual and gradient of the original LS problem
        r = b - A @ x
        grad = A.T @ r
        
        # If x is already an LS solution (or very close), retry with new random values
        if np.linalg.norm(grad) > 1e-6:
            break

    # We need to find a vector u that solves the system:
    # 1) A.T @ u + alpha * x = 0
    # 2) alpha = (r.T @ u - u.T @ u) / (x.T @ x)
    
    # The solution method depends on whether x is in the range of A.T
    # We can check this by solving a least-squares problem: min ||A.T @ z - x||^2
    z, residuals, _, _ = np.linalg.lstsq(A.T, x, rcond=None)
    
    if residuals.size > 0 and residuals[0] < 1e-9:
        # Case 1: x is in the range of A.T
        # We found two potential solutions for alpha, leading to two candidate u vectors.
        # u1 corresponds to alpha = 0.
        # u2 corresponds to the other root of a quadratic.
        u1 = np.zeros((m, 1))
        
        # Calculate the second candidate for alpha
        # From alpha * (||x||^2 + r.T@z + alpha*||z||^2) = 0
        if np.linalg.norm(z) > 1e-9:
            alpha2 = -(x.T @ x + r.T @ z) / (z.T @ z)
        else: # z is zero, so alpha can be anything. This case is unlikely for random data.
            alpha2 = 0
        u2 = -alpha2[0,0] * z

        # Choose the u that minimizes ||E||_F, which is equivalent to minimizing ||r-u||
        norm_r_minus_u1 = np.linalg.norm(r - u1)
        norm_r_minus_u2 = np.linalg.norm(r - u2)
        
        if norm_r_minus_u1 < norm_r_minus_u2:
            u = u1
        else:
            u = u2
            
    else:
        # Case 2: x is not in the range of A.T
        # This forces alpha = 0. The conditions become:
        # 1) A.T @ u = 0  (u is in the null space of A.T)
        # 2) r.T @ u - u.T @ u = 0
        # The optimal u that minimizes ||r-u|| under these constraints is the
        # orthogonal projection of r onto the null space of A.T.
        N = null_space(A.T)
        if N.shape[1] == 0: # Null space is empty
            u = np.zeros((m, 1))
        else:
            # Project r onto the null space of A.T
            u = N @ (N.T @ r)

    # Now that we have the optimal u, we can calculate E
    E = (r - u) @ x.T / (x.T @ x)

    # Calculate the rank of E
    rank_E = np.linalg.matrix_rank(E)

    # Print the "final equation" and the result
    print("The theoretical greatest possible rank of E is 1.")
    print("For a randomly generated example:")
    print(f"Dimensions: A is {m}x{n}, b is {m}x1, x is {n}x1.")
    print(f"The rank of the optimal perturbation matrix E is: {rank_E}")

    # You can uncomment the lines below to verify that x is indeed an LS solution for (A+E,b)
    # B = A + E
    # new_grad = B.T @ (B @ x - b)
    # print(f"Norm of new gradient ||(A+E).T @ ((A+E)x - b)||: {np.linalg.norm(new_grad)}")


if __name__ == '__main__':
    solve_greatest_rank_problem()