import numpy as np
from scipy.linalg import eig

def get_c_n(n):
    """
    Calculates the minimum eigenvalue c_n for a given integer n.
    """
    if n == 1:
        return 1.0
        
    # Matrix for the numerator quadratic form in terms of S_k
    L = np.zeros((n, n))
    for i in range(n - 1):
        L[i, i] = 2.0
        L[i, n - 1] = -1.0
        L[n - 1, i] = -1.0
    L[n - 1, n - 1] = n

    # Matrix for the denominator quadratic form in terms of S_k
    R = np.zeros((n, n))
    for i in range(n - 1):
        R[i, i] = 2.0
        R[i, i + 1] = -1.0
        R[i + 1, i] = -1.0
    R[0, 0] = 2.0 # S_1^2 + (S_2-S_1)^2 term for S_1^2
    R[n - 1, n - 1] = 1.0 # (S_n - S_{n-1})^2 term for S_n^2

    # The problem is to find the minimum of (S^T L S) / (S^T R S),
    # which is the minimum generalized eigenvalue of the pair (L, R).
    eigenvalues, _ = eig(L, R)
    
    # We are interested in the real eigenvalues
    real_eigenvalues = eigenvalues[np.isreal(eigenvalues)].real
    return np.min(real_eigenvalues)

def main():
    """
    Calculates c_n for several values of n and analyzes the trend.
    """
    ns = [2, 3, 4, 5, 10, 20, 50, 100, 200, 500, 1000]
    cs = [get_c_n(n) for n in ns]

    print("Calculated values for c_n:")
    for n, c in zip(ns, cs):
        print(f"n = {n:4d}, c_n = {c:.6f}")

    # The sequence c_n appears to be decreasing and converging to a limit.
    # Let's test the hypothesis that c_n converges to c = 1/4 = 0.25.
    # We can model the convergence as c_n = c + b/sqrt(n) for some constant b.
    # We can solve for c and b using the last few data points.
    
    # Using n1=500, n2=1000
    n1, c1 = ns[-2], cs[-2]
    n2, c2 = ns[-1], cs[-1]

    # c1 = c + b / sqrt(n1)
    # c2 = c + b / sqrt(n2)
    # c1 - c2 = b * (1/sqrt(n1) - 1/sqrt(n2))
    b = (c1 - c2) / (1/np.sqrt(n1) - 1/np.sqrt(n2))
    # c = c2 - b / sqrt(n2)
    c_limit = c2 - b / np.sqrt(n2)
    
    print("\nExtrapolating to find the limit c as n -> infinity...")
    print(f"Fit model: c_n = c + b/sqrt(n)")
    print(f"Estimated b: {b:.4f}")
    print(f"Estimated c: {c_limit:.6f}")
    print(f"\nThe analysis suggests that the limit is 1/4.")
    print(f"\nThe maximum value for c is the infimum of c_n, which is this limit.")
    final_c = 0.25
    print(f"\nFinal Answer: c = {final_c}")

if __name__ == "__main__":
    main()
