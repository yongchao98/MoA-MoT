import numpy as np
from scipy.optimize import minimize_scalar

def compute_and_show_rank():
    """
    This function constructs a specific least squares problem,
    calculates the minimal perturbation E, and prints its rank
    to demonstrate that a rank of 2 is achievable.
    """
    # 1. Define a problem instance (A, b, x) where x=e_1
    # We choose A to be 2x2, b to be 2x1.
    A = np.array([[1.0, 0.0], [1.0, 1.0]])
    b = np.array([-2.0, 0.0])
    x = np.array([1.0, 0.0])

    # Partition A for x=e_1 case
    a1 = A[:, 0]
    a2 = A[:, 1]
    
    # 2. Find the optimal vector v
    # The optimization problem can be reduced to a 1D minimization problem on a scalar s.
    # Let v = [v1, v2]^T. The constraint ||v||^2 + b^T*v = 0 implies v1^2+v2^2 - 2*v1 = 0.
    # We parameterize v using s=v1, so v2^2 = 2s-s^2 for s in [0, 2].
    
    # The objective function to minimize, as a function of s = v1.
    def objective_func(s):
        # Handle endpoints where v2 would be 0 or expression invalid.
        if s <= 0 or s >= 2:
            v2_squared = 0
        else:
            v2_squared = 2 * s - s**2
        
        v1 = s
        
        p = a1 - b
        term1_sq_norm = (v1 - p[0])**2 + (np.sqrt(v2_squared) - p[1])**2 if v2_squared>=0 else np.inf

        v_sq_norm = v1**2 + v2_squared
        if v_sq_norm == 0:
            return term1_sq_norm

        a2_T_v = a2[0] * v1 + a2[1] * np.sqrt(v2_squared) if v2_squared>=0 else 0
        term2_sq_norm = (a2_T_v**2) / v_sq_norm

        return term1_sq_norm + term2_sq_norm

    # Minimize the objective function to find the optimal s=v1
    res = minimize_scalar(objective_func, bounds=(0, 2), method='bounded')
    s_opt = res.x

    # Reconstruct the optimal vector v
    v1_opt = s_opt
    v2_opt = np.sqrt(2 * v1_opt - v1_opt**2)
    v = np.array([v1_opt, v2_opt])

    # 3. Calculate the perturbation matrix E
    # First column of E
    e1 = v + b - a1

    # Second column of E
    v_sq_norm = v.T @ v
    if v_sq_norm < 1e-9: # v is zero vector
        e2 = np.zeros_like(v)
    else:
        # For n=2, A2 is just a2. w is a scalar.
        w = -a2.T @ v
        e2 = v * w / v_sq_norm

    # Form the matrix E
    E = np.column_stack([e1, e2])

    # 4. Compute and print the rank of E
    rank_E = np.linalg.matrix_rank(E, tol=1e-8)
    
    print("For the constructed problem:")
    print(f"A = \n{A}")
    print(f"b = {b}")
    print(f"x = {x}")
    print(f"The computed optimal perturbation matrix E is:\n{E}")
    print(f"The rank of E is: {int(rank_E)}")

compute_and_show_rank()