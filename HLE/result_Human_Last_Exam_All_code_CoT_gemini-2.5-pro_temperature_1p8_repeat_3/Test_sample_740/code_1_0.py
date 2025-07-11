import numpy as np

def solve_for_lambda(n):
    """
    Calculates the values of lambda for which the given integral equation has no solution.

    Args:
        n (int): A positive integer, the exponent in the kernel of the integral equation.
    """
    if not isinstance(n, int) or n <= 0:
        print("Please provide a positive integer value for n.")
        return

    # A_ij = 1 / (i - j + n) for i, j from 0 to n-1
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            A[i, j] = 1.0 / (i - j + n)

    # Find the eigenvalues of matrix A.
    # The characteristic values lambda are the reciprocals of these eigenvalues.
    try:
        eigenvalues = np.linalg.eigvals(A)
        
        # Calculate lambda = 1 / mu for each eigenvalue mu
        # Check for non-zero eigenvalues to avoid division by zero.
        lambdas = [1/mu for mu in eigenvalues if abs(mu) > 1e-9]

        print(f"For n = {n}, the matrix A is:")
        print(A)
        print("\nThe eigenvalues (mu) of A are:")
        print(eigenvalues)
        print("\nThe values of lambda for which the equation has no solution are:")
        for l in lambdas:
            print(l)

    except np.linalg.LinAlgError:
        print(f"Could not compute eigenvalues for n = {n}.")


# You can set the integer value of 'n' here.
# For example, let's solve for n = 4.
n = 4
solve_for_lambda(n)
