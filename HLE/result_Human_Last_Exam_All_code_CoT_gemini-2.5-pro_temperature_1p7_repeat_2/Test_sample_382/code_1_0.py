import numpy as np
from scipy.optimize import root

def illustrate_e_rank():
    """
    This function illustrates that the rank of the minimum-norm matrix E is 1.
    It sets up a sample problem, solves for the perturbation E, and prints its rank.
    """
    # 1. Define a sample problem (A, b, x)
    # Let's choose m=3, n=2
    A = np.array([[1.0, 2.0],
                  [3.0, 4.0],
                  [5.0, 6.0]])
    b = np.array([[7.0],
                  [8.0],
                  [10.0]])
    # A non-zero vector x
    x = np.array([[1.0],
                  [1.0]])

    m, n = A.shape

    # Check if x is already a LS solution for A, b
    original_residual = b - A @ x
    original_normal_eq_val = A.T @ original_residual
    if np.linalg.norm(original_normal_eq_val) < 1e-9:
        print("x is already a least-squares solution. E is the zero matrix, rank is 0.")
        return

    # 2. Set up the non-linear system to find the vector u
    # The optimal E is of the form E = u @ x.T
    # The vector u is found by solving the system:
    # (||x||^2 * A.T - x @ r0.T) @ u + ||x||^2 * (u.T @ u) * x - A.T @ r0 = 0
    # where r0 = b - A @ x

    r0 = b - A @ x
    x_norm_sq = (np.linalg.norm(x)**2)
    M = x_norm_sq * A.T - x @ r0.T
    c = A.T @ r0

    # Define the function for the root finder
    # The function takes a flat array `u_flat` and returns the residuals of the system
    def func_for_u(u_flat):
        u = u_flat.reshape((m, 1))
        term1 = M @ u
        term2 = x_norm_sq * (u.T @ u) * x
        residual = term1 + term2 - c
        return residual.flatten()

    # 3. Solve for u using scipy.optimize.root
    # Initial guess for u
    u_initial_guess = np.zeros(m)
    solution = root(func_for_u, u_initial_guess, method='hybr')

    if not solution.success:
        print("Failed to find a solution for u.")
        return

    u = solution.x.reshape((m, 1))

    # 4. Construct E and compute its rank
    E = u @ x.T
    rank_E = np.linalg.matrix_rank(E)

    # 5. Print the results
    # As requested: "output each number in the final equation!"
    # The equation is E = u * x.T
    print("Given vector x:")
    print(x)
    print("\nSolved vector u:")
    print(u)
    print("\nResulting matrix E = u * x.T:")
    print(E)
    print(f"\nThe rank of the minimal perturbation matrix E is: {rank_E}")

if __name__ == '__main__':
    illustrate_e_rank()
