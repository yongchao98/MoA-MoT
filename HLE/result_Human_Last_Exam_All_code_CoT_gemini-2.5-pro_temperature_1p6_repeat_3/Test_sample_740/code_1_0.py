import numpy as np

def find_lambda_no_solution(n: int):
    """
    Calculates the values of λ for which the given integral equation has no solution.

    The integral equation u(x) = 1 + λ ∫[0,1] dy K(x,y) u(y) is converted into a
    system of linear equations (I - λM)c = b. The equation has no solutions
    for λ where 1/λ is an eigenvalue of the matrix M.

    Args:
        n: The integer exponent in the kernel K(x,y) = (xⁿ - yⁿ) / (x - y).

    Returns:
        A list of values for λ.
    """
    if n <= 0:
        print("Error: n must be a positive integer.")
        return []

    # The matrix M is derived from the kernel. Its elements are M_jk = 1/(n+j-k)
    # where j is the row index (0 to n-1) and k is the column index (0 to n-1).
    M = np.zeros((n, n))
    for j in range(n):
        for k in range(n):
            M[j, k] = 1.0 / (n + j - k)

    # The characteristic values λ are the reciprocals of the eigenvalues of M.
    try:
        eigenvalues = np.linalg.eigvals(M)
    except np.linalg.LinAlgError:
        print("Failed to compute eigenvalues for n =", n)
        return []

    # Filter out any zero eigenvalues as their reciprocal is undefined.
    # The matrix M is invertible, so it has no zero eigenvalues.
    lambda_values = 1.0 / eigenvalues

    # Print the results
    print(f"For n = {n}, the equation has no solutions for the following values of λ:")
    for val in lambda_values:
        # The eigenvalues of M can be complex for n>2, but since the integral operator
        # is self-adjoint, the eigenvalues of the operator (and its matrix representation)
        # must be real. Numerical precision might introduce small imaginary parts.
        # We print the real part.
        print(f"{val.real:.8f}")

    return lambda_values

# Example: Run the calculation for n = 3
if __name__ == '__main__':
    n_value = 3
    find_lambda_no_solution(n_value)

    # You can uncomment the line below to test with another value of n
    # find_lambda_no_solution(4)
