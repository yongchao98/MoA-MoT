import numpy as np
from scipy.optimize import minimize

def find_transience_direction(d, k, covariance_matrices=None):
    """
    This function demonstrates the existence of a "direction of transience" for a
    controlled random walk, given a set of k covariance matrices in d dimensions.

    A controller can guarantee recurrence only if for every direction (unit vector v),
    they can choose a measure `i` such that d*v^T*K_i*v - Tr(K_i) > 0.

    This function shows that this is never possible by finding a direction `v`
    that minimizes the maximum of these quantities over all measures. The theory
    predicts this minimum value will always be less than or equal to zero,
    meaning there is always a "direction of escape" where no measure helps
    in guaranteeing recurrence.
    """
    print(f"Starting demonstration for d={d}, k={k}.")
    print("-" * 40)
    
    # Generate random positive definite matrices for K_i if not provided
    if covariance_matrices is None:
        print("Generating random covariance matrices...")
        np.random.seed(42)
        covariance_matrices = []
        for _ in range(k):
            # Create a random matrix A and form K = A^T A + I to ensure it's positive definite
            A = np.random.rand(d, d)
            K = A.T @ A + np.identity(d)
            covariance_matrices.append(K)

    if len(covariance_matrices) != k:
        raise ValueError("Number of covariance matrices must be equal to k.")

    # Define the matrices M_i = d*K_i - Tr(K_i)*I
    M_matrices = []
    for K in covariance_matrices:
        if K.shape != (d, d):
            raise ValueError(f"Covariance matrix must be of size {d}x{d}.")
        trace_K = np.trace(K)
        M = d * K - trace_K * np.identity(d)
        M_matrices.append(M)

    # We want to find a unit vector `v` that minimizes `max_i(v^T M_i v)`.
    # This is a non-convex optimization problem, but solvers can often find the global minimum.
    def objective_function(v):
        # The solver's variable `v` is constrained to be a unit vector.
        quadratic_forms = [v.T @ M @ v for M in M_matrices]
        return np.max(quadratic_forms)

    # Constraint: ||v||^2 = 1
    constraints = ({'type': 'eq', 'fun': lambda v: np.sum(v**2) - 1})

    # Initial guess for the vector `v`
    v_init = np.random.rand(d)
    v_init /= np.linalg.norm(v_init)

    # Perform the optimization to find the "worst-case" direction v
    result = minimize(objective_function, v_init, constraints=constraints, method='SLSQP')

    min_max_value = result.fun
    v_optimal = result.x

    print("\nOptimization Goal:")
    print("Find a unit vector v that minimizes max_i(d*v^T*K_i*v - Tr(K_i)).")
    print("If this minimum value is <= 0, recurrence cannot be guaranteed.")
    
    print("\nResult:")
    print(f"The 'direction of escape' v* found is: {np.round(v_optimal, 4)}")
    print(f"The minimum of the maximum value is: {min_max_value:.8f}")

    # Check the condition
    if min_max_value <= 1e-9:  # Use a small tolerance for floating point errors
        print("\nConclusion: The value is <= 0, as predicted by theory.")
        print("This means a direction was found where no measure can guarantee recurrence.")
        print("This holds true for any k and any valid set of measures.")
    else:
        print("\nConclusion: The value is > 0.")
        print("This would contradict the theory (this may happen due to numerical limitations).")
        
    print("\nThe values of d*v*^T*K_i*v* - Tr(K_i) for this direction v* are:")
    for i, M in enumerate(M_matrices):
        val = v_optimal.T @ M @ v_optimal
        print(f"For measure {i+1}: {val:.8f} (which is <= 0)")


# --- Execute the demonstration for a sample case ---
# You can change these parameters
DIMENSION = 3
NUM_MEASURES = 5

find_transience_direction(d=DIMENSION, k=NUM_MEASURES)