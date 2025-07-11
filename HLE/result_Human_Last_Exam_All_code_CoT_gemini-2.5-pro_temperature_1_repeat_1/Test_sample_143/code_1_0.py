import numpy as np

def calculate_errors(w, X):
    """Calculates the squared error vector for a given w and data matrix X."""
    w = np.array(w) / np.linalg.norm(w) # Ensure w is a unit vector
    errors = []
    for x_i in X:
        x_i = np.array(x_i)
        # error_i_sq = ||x_i||^2 - (x_i . w)^2
        error_i_sq = np.dot(x_i, x_i) - np.dot(x_i, w)**2
        errors.append(error_i_sq)
    return np.array(errors)

def main():
    """
    Demonstrates the non-convexity of the Pareto front for d=3
    with non-negative data.
    """
    # Define the data matrix X (n=4, d=3) with non-negative entries
    X = np.array([
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.5, 0.5, 1.0/np.sqrt(2)]
    ])

    # Define two Pareto optimal solutions w_A and w_B
    w_A = np.array([1.0/np.sqrt(2), 1.0/np.sqrt(2), 0.0])
    w_B = np.array([0.0, 0.0, 1.0])

    # Calculate their corresponding error vectors
    e_A = calculate_errors(w_A, X)
    e_B = calculate_errors(w_B, X)

    # Calculate the midpoint of the error vectors
    e_mix = 0.5 * e_A + 0.5 * e_B

    # Define a third solution w_C which will dominate the midpoint
    w_C = np.array([0.5, 0.5, 1.0/np.sqrt(2)])
    e_C = calculate_errors(w_C, X)

    print("This script demonstrates that for d=3, the Pareto front can be non-convex even with non-negative data.")
    print("-" * 70)
    print("Data matrix X (all entries non-negative):")
    print(X)
    print("-" * 70)

    print(f"Error vector for w_A = {np.round(w_A, 3)}:")
    print(f"e_A = {np.round(e_A, 4)}")
    print(f"Equation: e_A = [||x1||^2-(x1.wA)^2, ||x2||^2-(x2.wA)^2, ||x3||^2-(x3.wA)^2, ||x4||^2-(x4.wA)^2]")
    print(f"e_A = [{1.0-0.5:.4f}, {1.0-0.5:.4f}, {1.0-0.0:.4f}, {1.0-0.5:.4f}] = [0.5000, 0.5000, 1.0000, 0.5000]")
    print("")

    print(f"Error vector for w_B = {np.round(w_B, 3)}:")
    print(f"e_B = {np.round(e_B, 4)}")
    print(f"Equation: e_B = [||x1||^2-(x1.wB)^2, ||x2||^2-(x2.wB)^2, ||x3||^2-(x3.wB)^2, ||x4||^2-(x4.wB)^2]")
    print(f"e_B = [{1.0-0.0:.4f}, {1.0-0.0:.4f}, {1.0-1.0:.4f}, {1.0-0.5:.4f}] = [1.0000, 1.0000, 0.0000, 0.5000]")
    print("")

    print("Convex combination (midpoint) of e_A and e_B:")
    print(f"e_mix = 0.5 * e_A + 0.5 * e_B")
    print(f"e_mix = {np.round(e_mix, 4)}")
    print(f"Equation: e_mix = 0.5 * [{e_A[0]:.4f}, {e_A[1]:.4f}, {e_A[2]:.4f}, {e_A[3]:.4f}] + 0.5 * [{e_B[0]:.4f}, {e_B[1]:.4f}, {e_B[2]:.4f}, {e_B[3]:.4f}]")
    print(f"e_mix = [{e_mix[0]:.4f}, {e_mix[1]:.4f}, {e_mix[2]:.4f}, {e_mix[3]:.4f}]")
    print("")

    print(f"Error vector for w_C = {np.round(w_C, 3)}:")
    print(f"e_C = {np.round(e_C, 4)}")
    print(f"Equation: e_C = [||x1||^2-(x1.wC)^2, ||x2||^2-(x2.wC)^2, ||x3||^2-(x3.wC)^2, ||x4||^2-(x4.wC)^2]")
    # ||w_C||^2 = 0.25+0.25+0.5 = 1, so it is a unit vector
    # x1.wC = 0.5, (x1.wC)^2 = 0.25. e1 = 1-0.25=0.75
    # x2.wC = 0.5, (x2.wC)^2 = 0.25. e2 = 1-0.25=0.75
    # x3.wC = 1/sqrt(2), (x3.wC)^2 = 0.5. e3 = 1-0.5=0.5
    # x4.wC = 0.5*0.5+0.5*0.5+(1/sqrt(2))*(1/sqrt(2)) = 0.25+0.25+0.5=1. (x4.wC)^2=1. e4=1-1=0
    print(f"e_C = [{1.0-0.25:.4f}, {1.0-0.25:.4f}, {1.0-0.5:.4f}, {1.0-1.0:.4f}] = [0.7500, 0.7500, 0.5000, 0.0000]")
    print("-" * 70)

    print("Comparing e_C with e_mix:")
    print(f"e_C   = {np.round(e_C, 4)}")
    print(f"e_mix = {np.round(e_mix, 4)}")
    print("\nConclusion:")
    print("e_C has errors equal to e_mix for the first three components, but a strictly smaller error for the fourth component:")
    print(f"e_C[0] = {e_C[0]:.4f} == e_mix[0] = {e_mix[0]:.4f}")
    print(f"e_C[1] = {e_C[1]:.4f} == e_mix[1] = {e_mix[1]:.4f}")
    print(f"e_C[2] = {e_C[2]:.4f} == e_mix[2] = {e_mix[2]:.4f}")
    print(f"e_C[3] = {e_C[3]:.4f} < e_mix[3] = {e_mix[3]:.4f}  <-- DOMINANCE")
    print("\nThis shows that e_C dominates e_mix. The existence of such a point proves that the Pareto front is non-convex.")
    print("Therefore, the condition of non-negative data is not sufficient for d=3.")
    print("The largest d for which the condition is sufficient is 2 (where convexity holds even without the condition).")

if __name__ == '__main__':
    main()
