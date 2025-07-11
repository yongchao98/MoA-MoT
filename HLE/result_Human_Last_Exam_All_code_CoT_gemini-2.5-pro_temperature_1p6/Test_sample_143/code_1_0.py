import numpy as np

def demonstrate_non_convexity():
    """
    Demonstrates that for d=3, the Pareto front can be non-convex even if all x_i >= 0.
    This implies scalarization is not sufficient.
    """
    
    # 1. We choose a counterexample for d=3 with n=2 non-negative vectors.
    x1 = np.array([1.0, 10.0, 0.0])
    x2 = np.array([10.0, 1.0, 0.0])
    
    print("--- Counterexample for d=3 ---")
    print(f"We use non-negative vectors x1 = {x1} and x2 = {x2}")
    print("We test if the set of achievable objectives K = { ((x1.w)^2, (x2.w)^2) | ||w||=1 } is convex.")
    
    # Define objective functions
    g = lambda w, x: np.dot(x, w)**2
    
    # 2. Find Point A on the Pareto front (by maximizing the first objective)
    w_A = x1 / np.linalg.norm(x1)
    A = np.array([g(w_A, x1), g(w_A, x2)])
    
    # 3. Find Point B on the Pareto front (by maximizing the second objective)
    w_B = x2 / np.linalg.norm(x2)
    B = np.array([g(w_B, x1), g(w_B, x2)])
    
    # 4. Find the midpoint M of the line segment AB
    M = (A + B) / 2
    
    # 5. Find another achievable point C
    w_C = np.array([1.0, 1.0, 0.0]) / np.sqrt(2)
    C = np.array([g(w_C, x1), g(w_C, x2)])

    print("\n--- Calculations ---")
    print("Point A (align w with x1):")
    print(f"g1 = (1 * {w_A[0]:.2f} + 10 * {w_A[1]:.2f})^2 = {A[0]:.2f}")
    print(f"g2 = (10 * {w_A[0]:.2f} + 1 * {w_A[1]:.2f})^2 = {A[1]:.2f}\n")
    
    print("Point B (align w with x2):")
    print(f"g1 = (1 * {w_B[0]:.2f} + 10 * {w_B[1]:.2f})^2 = {B[0]:.2f}")
    print(f"g2 = (10 * {w_B[0]:.2f} + 1 * {w_B[1]:.2f})^2 = {B[1]:.2f}\n")
    
    print("Midpoint M of AB:")
    print(f"g1 = ({A[0]:.2f} + {B[0]:.2f}) / 2 = {M[0]:.2f}")
    print(f"g2 = ({A[1]:.2f} + {B[1]:.2f}) / 2 = {M[1]:.2f}\n")
    
    print("Point C (w proportional to (1,1,0)):")
    print(f"g1 = (1 * {w_C[0]:.2f} + 10 * {w_C[1]:.2f})^2 = {C[0]:.2f}")
    print(f"g2 = (10 * {w_C[0]:.2f} + 1 * {w_C[1]:.2f})^2 = {C[1]:.2f}\n")

    # 6. Conclusion from comparison
    print("--- Conclusion ---")
    print(f"We compare C({C[0]:.2f}, {C[1]:.2f}) with M({M[0]:.2f}, {M[1]:.2f}).")
    if np.all(C > M):
        print("Point C dominates M. This shows the achievable objective set is not convex.")
        print("Because the front is not convex, scalarization is not guaranteed to find all Pareto optimal solutions.")
        print("Since a counterexample exists for d=3, but d=2 is known to be sufficient, the largest d is 2.")
    else:
        print("The test was inconclusive.")

demonstrate_non_convexity()