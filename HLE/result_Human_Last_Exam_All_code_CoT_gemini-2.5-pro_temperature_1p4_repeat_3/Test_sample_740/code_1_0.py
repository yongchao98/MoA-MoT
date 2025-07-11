import numpy as np

def find_lambdas_for_no_solution(n):
    """
    Finds the values of λ for which the given integral equation has no solution, for a specified integer n.
    """
    if not isinstance(n, int) or n < 1:
        print("Error: n must be a positive integer.")
        return []

    # Construct the matrix A, where A[j,k] = 1/(j+n-k)
    A = np.fromfunction(lambda j, k: 1.0 / (j + n - k), (n, n), dtype=int)

    # The integral operator is self-adjoint with respect to the L2 inner product.
    # Its eigenvalues are real. The matrix A, which represents the operator in a
    # non-orthogonal basis, is not symmetric but can be shown to have real eigenvalues.
    try:
        eigenvalues, eigenvectors = np.linalg.eig(A)
    except np.linalg.LinAlgError:
        print(f"Error: Eigenvalue computation failed for n={n}.")
        return []

    no_solution_lambdas = []
    
    for i in range(n):
        mu = eigenvalues[i]
        
        # The eigenvector c is the i-th column of the eigenvectors matrix
        c_vec = eigenvectors[:, i]
        
        # The condition for no solution is that the integral of the corresponding
        # eigenfunction v(x) is non-zero. For an eigenfunction proportional to
        # Σ c_k * x^(n-1-k), the integral ∫v(x)dx is proportional to Σ c_k / (n-k).
        integral_check = np.sum(c_vec / np.arange(n, 0, -1))
            
        # We check if the integral condition is met (i.e., not close to zero)
        if abs(integral_check) > 1e-9:
            # If the condition holds, λ = 1/μ is a value for which no solution exists.
            # An eigenvalue μ=0 is highly unlikely for this type of matrix. If it occurs,
            # λ would be infinite, which is not a valid characteristic value.
            if abs(mu) > 1e-9:
                lambda_val = 1.0 / mu
                no_solution_lambdas.append(lambda_val)
                
    # Since eigenvalues of a real matrix can appear as complex conjugate pairs,
    # and we expect real values, we take the real part. Any imaginary part
    # should be due to numerical floating-point inaccuracies.
    return sorted([l.real for l in no_solution_lambdas])

# --- Execution ---

# Set the value of n for the equation u(x) = 1 + λ ∫₀¹ dy (xⁿ - yⁿ) / (x - y) u(y)
# You can change this value to solve for different n.
n = 4

print(f"This script calculates the values of λ for which the integral equation has no solutions.")
print(f"The values of λ depend on the integer 'n'.")
print(f"Running for n = {n}\n")

result_lambdas = find_lambdas_for_no_solution(n)

if result_lambdas:
    print(f"For n = {n}, the equation has no solution for the following values of λ:")
    # Print the numbers in the final result, rounded for clarity
    print(np.round(result_lambdas, 6))
else:
    print(f"Could not determine the values of λ for n = {n}.")
    
# For example, to find the values for n=1 and n=2:
print("\n--- Examples ---")
print("For n = 1, λ = ", np.round(find_lambdas_for_no_solution(1), 6))
print("For n = 2, λ = ", np.round(find_lambdas_for_no_solution(2), 6))
