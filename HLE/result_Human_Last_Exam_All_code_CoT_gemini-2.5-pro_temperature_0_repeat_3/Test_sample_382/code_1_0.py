import numpy as np
from scipy.optimize import minimize

def solve_and_find_rank():
    """
    This function sets up a least-squares problem, finds the minimum-norm
    perturbation E such that a given vector x is a solution, and computes
    the rank of E.
    """
    # 1. Define the problem matrices and vectors
    # We choose A, b, x such that x is NOT a least-squares solution for Ax=b.
    # A is a 3x2 matrix, b is a 3x1 vector, x is a 2x1 vector.
    m, n = 3, 2
    A = np.array([[1.0, 2.0], 
                  [3.0, 4.0], 
                  [5.0, 6.0]])
    
    b = np.array([1.0, 1.0, 1.0])
    
    x = np.array([1.0, 1.0])

    # Check if x is a LS solution for the original problem
    original_residual = A @ x - b
    original_normal_eq_val = A.T @ original_residual
    
    print("--- Problem Setup ---")
    print("Matrix A:\n", A)
    print("Vector b:\n", b)
    print("Vector x:\n", x)
    print("\nChecking if x is already a LS solution for (A,b):")
    print("A.T @ (A@x - b) =", original_normal_eq_val)
    if np.allclose(original_normal_eq_val, 0):
        print("x is already a LS solution. E should be the zero matrix (rank 0).")
    else:
        print("x is NOT a LS solution. E should be a non-zero, rank-1 matrix.")

    # 2. Define the optimization functions
    # The variable for the optimizer is the flattened matrix E.
    
    # Objective function: minimize the Frobenius norm of E, which is ||E||_F^2
    def objective(E_vec):
        return np.sum(E_vec**2)

    # Constraint function: the normal equations for the perturbed system must be zero.
    # (A+E)^T * ((A+E)x - b) = 0
    def constraint(E_vec):
        E = E_vec.reshape((m, n))
        B = A + E
        residual = B @ x - b
        # The constraint is that this vector must be the zero vector
        return B.T @ residual

    # 3. Solve the constrained optimization problem
    # Initial guess for E is the zero matrix
    E0_vec = np.zeros(m * n)

    # Define the constraint for the solver
    cons = {'type': 'eq', 'fun': constraint}

    # Use the SLSQP (Sequential Least Squares Programming) method
    result = minimize(objective, E0_vec, method='SLSQP', constraints=cons, tol=1e-9)

    if not result.success:
        print("\nOptimization failed:", result.message)
        return

    # 4. Analyze the result
    # Reshape the solution vector back into a matrix
    E_sol = result.x.reshape((m, n))

    # Calculate the rank of the resulting matrix E
    # Use a small tolerance to account for numerical inaccuracies
    rank = np.linalg.matrix_rank(E_sol, tol=1e-7)

    print("\n--- Solution ---")
    print("The optimizer found the minimum-norm matrix E:\n", E_sol)
    print("\nFrobenius norm of E:", np.linalg.norm(E_sol))
    print("\nRank of the resulting matrix E:", rank)
    print("\nThe greatest possible rank for E is 1.")


if __name__ == '__main__':
    solve_and_find_rank()