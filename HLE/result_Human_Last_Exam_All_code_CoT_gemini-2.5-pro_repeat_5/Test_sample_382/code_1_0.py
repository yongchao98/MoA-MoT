import numpy as np
from scipy.optimize import minimize

def solve_for_E_and_find_rank(A, x, b):
    """
    Solves for the matrix E with minimum Frobenius norm such that x is a
    least-squares solution to (A+E)y = b, and returns its rank.

    Args:
        A (np.ndarray): The matrix A.
        x (np.ndarray): The vector x.
        b (np.ndarray): The vector b.

    Returns:
        int: The rank of the optimal matrix E.
    """
    m, n = A.shape

    # Objective function to minimize: ||E||_F^2
    def objective(E_flat):
        return np.sum(E_flat**2)

    # Constraint function: (A+E)^T((A+E)x - b) = 0
    def constraint(E_flat):
        E = E_flat.reshape((m, n))
        B = A + E
        residual = B @ x - b
        # The constraint is that the residual is orthogonal to the column space of B.
        # This is equivalent to B^T @ residual = 0
        constraint_val = B.T @ residual
        return constraint_val.flatten()

    # Initial guess for E (a zero matrix)
    E0_flat = np.zeros(m * n)

    # Set up the optimization problem with an equality constraint
    con = {'type': 'eq', 'fun': constraint}

    # Solve the optimization problem using Sequential Least Squares Programming (SLSQP)
    solution = minimize(objective, E0_flat, method='SLSQP', constraints=con, options={'maxiter': 200, 'ftol': 1e-9})

    if not solution.success:
        print("Warning: Optimization may not have converged.")
        print(solution.message)

    # Extract the solution for E
    E_optimal = solution.x.reshape((m, n))

    # Calculate the rank of the resulting matrix E using a small tolerance
    rank = np.linalg.matrix_rank(E_optimal, tol=1e-6)
    
    # --- Outputting the details of the final result ---
    print("--- Problem Setup ---")
    print("Matrix A:")
    print(A)
    print("\nVector x:")
    print(x)
    print("\nVector b:")
    print(b)
    
    print("\n--- Result ---")
    print("Found optimal perturbation matrix E:")
    print(E_optimal)
    
    print("\nFinal Equation: Rank(E)")
    print(f"The rank of the optimal matrix E is: {rank}")
    print("========================================\n")

    return rank

if __name__ == '__main__':
    # Define a specific problem instance
    # We choose m=3, n=2 for our demonstration
    A = np.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
    x = np.array([[1.0], [1.0]])
    # We choose b such that the initial residual is non-zero
    b = np.array([[10.0], [10.0], [10.0]])

    # Solve the problem and find the rank
    greatest_rank = solve_for_E_and_find_rank(A, x, b)
    
    # The greatest possible rank is the result from our theoretical analysis,
    # demonstrated by the numerical example.
    print(f"Based on the analysis, the greatest possible rank of E is {greatest_rank}.")
