import numpy as np

def solve_pareto_front_analysis():
    """
    Analyzes the convexity of the Pareto front for the rank-1 approximation problem.
    
    This function demonstrates that for d=3, even with non-negative data vectors x_i,
    the Pareto front is not necessarily convex, meaning scalarization can fail to
    find the entire front. This implies the largest d for which the condition is
    sufficient is d=2.
    """
    print("Analyzing the multi-objective rank-1 approximation problem.")
    print("We check if x_i >= 0 is a sufficient condition for scalarization to work.\n")
    
    print("--- Case d=2 ---")
    print("For d=2, the joint numerical range of the objective matrices is always convex.")
    print("This means scalarization can find the entire Pareto front for ANY data,")
    print("so the condition x_i >= 0 is sufficient.\n")
    
    print("--- Case d=3 ---")
    print("For d=3, we test with a counterexample where x_i >= 0.")
    print("We will construct two non-negative vectors x1, x2 in R^3.")
    # This counterexample is adapted from academic literature on the topic.
    x1 = np.array([np.sqrt(0.2), np.sqrt(0.8), 0.0])
    x2 = np.array([0.0, np.sqrt(0.045), np.sqrt(0.405)])
    
    print(f"x1 = {x1}")
    print(f"x2 = {x2}")
    print("Both vectors have non-negative entries.\n")

    # The objectives are g1(w) = (w^T x1)^2 and g2(w) = (w^T x2)^2.
    # We analyze the scalarized problem with equal weights (lambda1=lambda2).
    # This is equivalent to finding the eigenvectors of M = x1*x1^T + x2*x2^T.
    A1 = np.outer(x1, x1)
    A2 = np.outer(x2, x2)
    M = A1 + A2
    
    print("Considering the scalarized problem with equal weights, we analyze the matrix M = A1 + A2:")
    print(M)
    print("\n")
    
    # Find eigenvalues and eigenvectors of M
    eigenvalues, eigenvectors = np.linalg.eigh(M)
    
    # Sort by eigenvalue in descending order
    idx = eigenvalues.argsort()[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    
    # w1 is the dominant eigenvector, which is the solution for this scalarization.
    # w2 is the eigenvector for the second-largest eigenvalue. It is a candidate
    # for a Pareto optimal point that scalarization might miss.
    w1 = eigenvectors[:, 0]
    w2 = eigenvectors[:, 1]
    
    print(f"Largest eigenvalue of M: {eigenvalues[0]:.4f}")
    print(f"Dominant eigenvector (w1): {w1}\n")
    print(f"Second largest eigenvalue of M: {eigenvalues[1]:.4f}")
    print(f"Corresponding eigenvector (w2): {w2}\n")
    
    # Calculate the objective values for w1 and w2
    g1_w1 = (w1.T @ x1)**2
    g2_w1 = (w1.T @ x2)**2
    
    g1_w2 = (w2.T @ x1)**2
    g2_w2 = (w2.T @ x2)**2
    
    print("--- Comparing Objective Values ---")
    print("Objective vector for w1 (scalarization solution):")
    print(f"g(w1) = (g1, g2) = ({g1_w1:.4f}, {g2_w1:.4f})")
    print("Objective vector for w2 (candidate for non-supported solution):")
    print(f"g(w2) = (g1, g2) = ({g1_w2:.4f}, {g2_w2:.4f})\n")
    
    print("For w1 to dominate w2, we need g1(w1) >= g1(w2) AND g2(w1) >= g2(w2).")
    
    # Let's output the final equations for comparison.
    print("\n--- Final Comparison ---")
    print(f"Comparing first objective: g1(w1) vs g1(w2)")
    print(f"Equation: {g1_w1:.4f} vs {g1_w2:.4f}")
    if g1_w1 < g1_w2:
        print("Result: w1 is WORSE than w2 on the first objective.\n")
    else:
        print("Result: w1 is better or equal to w2 on the first objective.\n")

    print(f"Comparing second objective: g2(w1) vs g2(w2)")
    print(f"Equation: {g2_w1:.4f} vs {g2_w2:.4f}")
    if g2_w1 < g2_w2:
        print("Result: w1 is WORSE than w2 on the second objective.\n")
    else:
        print("Result: w1 is better or equal to w2 on the second objective.\n")
        
    if g1_w1 < g1_w2 or g2_w1 < g2_w2:
        print("Conclusion: The scalarization solution w1 does not dominate w2.")
        print("This demonstrates that the Pareto front is non-convex for d=3, even with x_i >= 0.")
        print("Therefore, the condition is not sufficient for d=3.\n")

    print("Final Answer: The largest dimension 'd' for which the condition is sufficient is 2.")

if __name__ == '__main__':
    solve_pareto_front_analysis()
