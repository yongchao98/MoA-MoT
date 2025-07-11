import numpy as np
import cvxpy as cp
from itertools import product

def solve_for_z():
    """
    Calculates the smallest z for a specific 5x5 correlation matrix A,
    which serves as a lower bound for the Grothendieck constant K_G.
    """
    # 1. Define the dimension n
    n = 5

    # 2. Define a specific correlation matrix A (from the pentagon graph)
    # This matrix is known to be a challenging instance for this type of problem.
    c1 = np.cos(2 * np.pi / 5)
    c2 = np.cos(4 * np.pi / 5)
    A = np.array([
        [1, c1, c2, c2, c1],
        [c1, 1, c1, c2, c2],
        [c2, c1, 1, c1, c2],
        [c2, c2, c1, 1, c1],
        [c1, c2, c2, c1, 1]
    ])

    print("The problem is to find the smallest 'z' for a given A.")
    print("We will solve it for a specific 5x5 matrix A, which gives a lower bound for the true answer K_G.\n")

    # 3. Generate the basis matrices for "nice" matrices
    # These are matrices of the form s*s^T, where s has entries in {-1, 1}.
    # We only need one from each pair (s, -s), so there are 2^(n-1) unique matrices.
    s_vectors = list(product([-1, 1], repeat=n))
    unique_s_tuples = {tuple(s) for s in s_vectors}
    # Halve the set by removing symmetric vectors
    seen_s_tuples = set()
    unique_s_vectors = []
    for s_tuple in unique_s_tuples:
        neg_s_tuple = tuple([-x for x in s_tuple])
        if neg_s_tuple not in seen_s_tuples:
            unique_s_vectors.append(np.array(s_tuple))
            seen_s_tuples.add(s_tuple)
    
    basis_matrices = [np.outer(s, s) for s in unique_s_vectors]
    num_basis = len(basis_matrices)

    # 4. Set up and solve the Semidefinite Program (SDP)
    # Variables: z (scalar) and lambdas (weights for the convex combination)
    z = cp.Variable()
    lambdas = cp.Variable(num_basis, nonneg=True)

    # The "nice" matrix B is a convex combination of the basis matrices
    B = cp.sum([lambdas[k] * basis_matrices[k] for k in range(num_basis)])

    # Constraints for the optimization
    constraints = [
        cp.sum(lambdas) == 1,   # The weights must sum to 1 (convex combination)
        z * B - A >> 0          # The matrix z*B - A must be positive semidefinite
    ]

    # The objective is to find the minimum possible z
    objective = cp.Minimize(z)

    # Define and solve the problem
    problem = cp.Problem(objective, constraints)
    problem.solve()

    # 5. Output the results and the final equation numbers
    if problem.status == cp.OPTIMAL:
        z_val = problem.value
        lambdas_val = lambdas.value

        # Construct the optimal matrices B and C
        B_val = np.zeros((n, n))
        for k in range(num_basis):
            B_val += lambdas_val[k] * basis_matrices[k]

        C_val = z_val * B_val - A
        
        print(f"The smallest 'z' for this specific matrix A is: {z_val:.4f}\n")
        print("This calculated 'z' is a lower bound for the universal constant K_G.")
        print("-" * 50)
        print("The final equation is A = z*B - C, where:\n")
        
        print("z =")
        print(f"{z_val:.4f}\n")
        
        print("A (the given correlation matrix) =")
        print(np.round(A, 4))
        print("\n")

        print("The term z*B in the equation =")
        print(np.round(z_val * B_val, 4))
        print("\n")

        print("C = z*B - A (must be positive semidefinite) =")
        print(np.round(C_val, 4))
        print("\n")

        c_eigenvalues = np.linalg.eigvalsh(C_val)
        print("Verification: Eigenvalues of C (should be non-negative):")
        # Eigenvalues might be slightly negative due to solver precision, e.g., -1e-9
        print(np.round(c_eigenvalues, 8))

    else:
        print("The optimization problem could not be solved.")

if __name__ == '__main__':
    solve_for_z()