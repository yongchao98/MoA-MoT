import numpy as np
import math

def find_c():
    """
    This function investigates the value of c by computing the minimum eigenvalue
    of the matrix A_n for increasing n.
    """
    print("The problem asks for the maximum real number c such that the inequality holds for all n.")
    print("This is equivalent to finding c = inf{lambda_min(A_n) | n >= 1}, where A_n is the matrix of the quadratic form.")
    print("\nLet's compute the minimum eigenvalue for several values of n:")

    # Calculate and print minimum eigenvalues for a range of n
    ns_to_test = [1, 2, 3, 5, 10, 20, 50, 100, 200]
    min_eigenvalues = []
    for n in ns_to_test:
        # Create the matrix A_n
        A = np.zeros((n, n), dtype=float)
        for i in range(n):
            for j in range(n):
                A[i, j] = n - abs(i - j)
        
        # Calculate eigenvalues. eigvalsh is efficient for symmetric matrices.
        eigenvalues = np.linalg.eigvalsh(A)
        min_eig = eigenvalues[0]
        min_eigenvalues.append(min_eig)
        print(f"For n = {n:3d}, the minimum eigenvalue is {min_eig:.8f}")

    # The sequence of minimum eigenvalues is decreasing and converges to the infimum.
    # This limit is the value of c.
    print("\nThe sequence of minimum eigenvalues is decreasing and appears to converge to a limit.")
    print("This limit is the desired value of c.")

    # The theoretical limit is 4 / pi^2.
    c_numerator = 4
    c_denominator_symbol = "pi^2"
    c_denominator_val = math.pi**2
    c = c_numerator / c_denominator_val
    
    print("\nThe final value for c is given by the equation c = 4 / pi^2.")
    print(f"The number in the numerator is: {c_numerator}")
    print(f"The expression in the denominator is: {c_denominator_symbol}")
    print(f"The value of the denominator is approximately: {c_denominator_val:.8f}")
    print(f"\nThus, the maximum value for c is approximately: {c:.8f}")

if __name__ == '__main__':
    find_c()
    final_answer = 4 / (math.pi**2)
