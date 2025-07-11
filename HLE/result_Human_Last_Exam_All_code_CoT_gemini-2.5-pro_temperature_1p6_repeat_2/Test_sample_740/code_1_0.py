import numpy as np

def find_lambda_for_no_solution(n):
    """
    Calculates the values of λ for which the given integral equation has no solution.

    Args:
        n (int): The integer exponent in the kernel of the integral equation.

    Returns:
        list: A list of values of λ.
    """
    # According to the derivation, the values of λ are the reciprocals of the
    # eigenvalues of the matrix M, where M_kj = 1/(n+k-j) for k,j = 0..n-1.

    # Step 1: Construct the matrix M.
    M = np.zeros((n, n))
    for k in range(n):
        for j in range(n):
            M[k, j] = 1.0 / (n + k - j)

    # Step 2: Find the eigenvalues of M.
    try:
        eigenvalues = np.linalg.eigvals(M)
    except np.linalg.LinAlgError:
        print(f"Eigenvalue computation failed for n={n}.")
        return []

    # Step 3: λ are the reciprocals of the eigenvalues.
    # We filter out eigenvalues close to zero to avoid division by zero errors.
    lambda_values = []
    for mu in eigenvalues:
        if not np.isclose(mu, 0):
            lambda_values.append(1.0 / mu)

    return lambda_values

def main():
    """
    Main function to execute the solution.
    The problem does not specify the integer n. We demonstrate the solution for n=2.
    """
    n = 2
    print(f"For n = {n}, the values of λ for which the equation has no solution are:")

    # For n=2, we can also solve analytically.
    # The matrix M is [[1/2, 1], [1/3, 1/2]].
    # The eigenvalues μ satisfy (1/2-μ)^2 - 1/3 = 0, so μ = 1/2 ± 1/√3.
    # The values for λ are λ = 1/μ = -6 ± 4√3.
    # The numerical values are approximately 0.928 and -12.928.
    
    # We will now compute them using our function.
    lambda_vals = find_lambda_for_no_solution(n)
    
    # Sort the values for consistent output
    lambda_vals.sort()
    
    for val in lambda_vals:
        print(f"{val:.8f}")

if __name__ == '__main__':
    main()
