import numpy as np

def find_lambda_for_no_solution(n):
    """
    Calculates the values of lambda for which the integral equation has no solutions.

    Args:
        n (int): The integer parameter 'n' in the integral equation. Must be >= 1.
    """
    if not isinstance(n, int) or n < 1:
        print("Error: n must be a positive integer.")
        return

    # Create the n x n matrix A
    # A_jk = 1 / (j + n - k) where j, k are 0-indexed.
    A = np.zeros((n, n))
    for j in range(n):
        for k in range(n):
            A[j, k] = 1.0 / (j + n - k)

    # Find the eigenvalues of matrix A.
    # The matrix A may not be symmetric, so we use eigvals for general matrices.
    try:
        eigenvalues = np.linalg.eigvals(A)
    except np.linalg.LinAlgError:
        print("Could not compute the eigenvalues for the matrix A.")
        return

    # The values of lambda are the reciprocals of the eigenvalues of A.
    # Eigenvalues of A are guaranteed to be non-zero since A is invertible.
    lambda_values = 1.0 / eigenvalues

    print(f"For n = {n}, the integral equation has no solution for the following values of lambda:")
    for val in lambda_values:
        # NumPy might return complex numbers with a tiny imaginary part for real results
        if np.iscomplex(val):
            print(f"{val.real: .8f} + {val.imag: .8f}j")
        else:
            print(f"{val.real: .8f}")

if __name__ == '__main__':
    try:
        # Get the value of n from the user
        n_input = int(input("Enter the integer value for n: "))
        find_lambda_for_no_solution(n_input)
    except ValueError:
        print("Invalid input. Please enter an integer.")
