import numpy as np
import sys

def solve_for_lambda(n):
    """
    Calculates the values of lambda for which the integral equation has no solution
    for a given integer n >= 1.
    """
    if not isinstance(n, int) or n < 1:
        print("Error: n must be a positive integer.", file=sys.stderr)
        return

    # Construct the matrix A
    A = np.zeros((n, n))
    for k in range(n):
        for j in range(n):
            A[k, j] = 1.0 / (k + n - j)

    # Calculate the eigenvalues of A
    try:
        eigenvalues = np.linalg.eigvals(A)
    except np.linalg.LinAlgError:
        print(f"Error: Eigenvalue computation failed for n={n}.", file=sys.stderr)
        return

    # Lambda values are the reciprocals of the eigenvalues
    # Filter out potential division by zero if an eigenvalue is 0
    lambdas = [1/mu for mu in eigenvalues if abs(mu) > 1e-9]

    print(f"For n = {n}, the equation has no solution for the following values of lambda:")
    for val in lambdas:
        # The eigenvalues of this real, non-symmetric matrix can be complex
        if np.iscomplex(val):
            print(f"{val.real:.8f} + {val.imag:.8f}j")
        else:
            print(f"{val.real:.8f}")

if __name__ == '__main__':
    # You can change the value of n here.
    # For example, n=1, n=2, n=3, etc.
    try:
        n_val = int(input("Enter the value of n (a positive integer): "))
        solve_for_lambda(n_val)
    except (ValueError, TypeError):
        print("Invalid input. Please enter a positive integer for n.", file=sys.stderr)
