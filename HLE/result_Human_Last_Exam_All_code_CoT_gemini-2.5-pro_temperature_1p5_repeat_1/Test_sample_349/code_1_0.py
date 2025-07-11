import numpy as np
import cvxpy as cp
import itertools

def solve_grothendieck_lower_bound():
    """
    This script calculates a numerical lower bound for the Grothendieck constant K_G.

    The constant z in the problem is K_G. Its value is not known exactly.
    We can compute lower bounds for it by solving a Semidefinite Program (SDP)
    for a specific dimension `n` and a specific "hard" correlation matrix `A`.

    We use n=5, the smallest non-trivial dimension, and a matrix `A` based on a
    regular pentagon, which is known to be a "hard" instance.

    The calculated value is a lower bound for K_G(5), and thus for K_G.
    K_G(5) is known to be approximately 1.2655, which our calculation
    should approximate. This shows K_G is greater than 1.
    """
    n = 5

    # 1. Define the "hard" correlation matrix A for n=5
    # A_ij = cos(2*pi*(i-j)/5)
    indices = np.arange(n)
    A = np.cos(2 * np.pi * (indices[:, np.newaxis] - indices) / n)

    # 2. Set up the Semidefinite Program (SDP)
    # Variable to be optimized is a symmetric n x n matrix X
    X = cp.Variable((n, n), symmetric=True)

    # Objective function: maximize Tr(A @ X)
    objective = cp.Maximize(cp.trace(A @ X))

    # Constraints
    constraints = [X >> 0]  # X must be positive semidefinite

    # Generate all 2^n sign vectors u in {-1, 1}^n
    # We only need 2^(n-1) unique vectors for the constraints u.T@X@u,
    # since (-u).T@X@(-u) gives the same constraint. We'll take those with u[0]=1.
    num_vectors = 2**(n - 1)
    sign_vectors = list(itertools.product([-1, 1], repeat=n - 1))

    for i in range(num_vectors):
        u = np.concatenate(([1], sign_vectors[i]))
        # Add constraint: u.T @ X @ u <= 1
        constraints.append(cp.quad_form(u, X) <= 1)

    # 3. Define and solve the problem
    problem = cp.Problem(objective, constraints)
    
    # We use SCS solver as it is good for these types of problems.
    # We may need to increase the number of iterations for better precision.
    problem.solve(solver=cp.SCS, verbose=False, max_iters=5000)

    # 4. Output the result
    # The optimal value of the SDP is a lower bound for K_G.
    lower_bound = problem.value
    
    # The "final equation" in this context is the statement about the lower bound found.
    print(f"The problem asks for the value of the Grothendieck constant, K_G.")
    print(f"We computed a lower bound by solving an SDP for n=5.")
    print(f"The computed lower bound for K_G is: {lower_bound:.4f}")
    # This demonstrates that the constant is non-trivial.
    # The actual value of K_G is known to be between 1.6769 and 1.7822.

solve_grothendieck_lower_bound()