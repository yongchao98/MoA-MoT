import numpy as np

def demonstrate_non_convexity():
    """
    Demonstrates that the Pareto front for the rank-1 approximation problem
     can be non-convex even if all data points x_i have non-negative entries.
    This implies that the weighted sum (scalarization) method cannot find all
    Pareto optimal solutions.
    """
    # 1. Define the data matrix X with non-negative entries
    # Let's choose d=2, n=2
    epsilon = 0.5
    X = np.array([[1.0, epsilon], [epsilon, 1.0]])
    x1 = X[0, :]
    x2 = X[1, :]
    
    print(f"Using counterexample with X =\n{X}\n")

    # 2. Define a function to calculate the objective values
    def calculate_objectives(w, X):
        """Calculates f_i(w) = (x_i^T w)^2 for i=1, ..., n"""
        objectives = []
        for i in range(X.shape[0]):
            xi = X[i, :]
            obj_val = (xi @ w)**2
            objectives.append(obj_val)
        return np.array(objectives)

    # 3. Find two points P_A and P_B on the Pareto front
    
    # P_A corresponds to the vector w_A that maximizes the first objective f_1
    # This vector is aligned with x1
    w_A = x1 / np.linalg.norm(x1)
    P_A = calculate_objectives(w_A, X)

    # P_B corresponds to the vector w_B that maximizes the second objective f_2
    # This vector is aligned with x2
    w_B = x2 / np.linalg.norm(x2)
    P_B = calculate_objectives(w_B, X)

    print(f"Point P_A on Pareto front: {P_A}")
    print(f"Point P_B on Pareto front: {P_B}")

    # 4. Calculate the midpoint of the chord connecting P_A and P_B
    P_M = (P_A + P_B) / 2.0
    print(f"Midpoint P_M of the chord P_A-P_B: {P_M}")

    # 5. Find a third point P_C that might dominate P_M
    # Let's choose w_C to be a direction "in between" w_A and w_B
    w_C = np.array([1.0, 1.0]) / np.sqrt(2)
    P_C = calculate_objectives(w_C, X)
    print(f"A third achievable point P_C: {P_C}")

    # 6. Check if P_C dominates P_M
    # This demonstrates that the Pareto front is non-convex
    print("\n--- Checking for Dominance ---")
    is_dominated = (P_C[0] > P_M[0]) and (P_C[1] > P_M[1])
    
    print("For the front to be non-convex, P_C must dominate P_M.")
    print(f"Is P_C's 1st objective > P_M's 1st objective? \
({P_C[0]:.4f} > {P_M[0]:.4f})")
    print(f"Final equation check: {P_C[0]} > {P_M[0]} is {P_C[0] > P_M[0]}")
    
    print(f"Is P_C's 2nd objective > P_M's 2nd objective? \
({P_C[1]:.4f} > {P_M[1]:.4f})")
    print(f"Final equation check: {P_C[1]} > {P_M[1]} is {P_C[1] > P_M[1]}")

    if is_dominated:
        print("\nConclusion: P_C dominates P_M. The Pareto front is non-convex.")
        print("Therefore, the condition x_i >= 0 is not sufficient for d=2.")
    else:
        print("\nConclusion: The chosen P_C does not dominate P_M.")

demonstrate_non_convexity()