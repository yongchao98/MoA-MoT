import numpy as np

def solve_for_lambda(n):
    """
    Finds the values of lambda for which the integral equation has no solutions,
    for a given integer n.
    
    The equation is: u(x) = 1 + lambda * integral from 0 to 1 of K(x,y)u(y) dy
    with K(x,y) = (x^n - y^n) / (x - y).
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return

    # 1. Construct the matrix A
    # A_ji = 1 / (n + j - i) for j, i in {0, ..., n-1}
    A = np.zeros((n, n))
    for j in range(n):
        for i in range(n):
            A[j, i] = 1.0 / (n + j - i)

    # 2. Find eigenvalues and eigenvectors of A
    # The operator is self-adjoint, so eigenvalues must be real.
    # We take the real part to handle potential small numerical imaginary parts.
    eigenvalues, eigenvectors = np.linalg.eig(A)
    eigenvalues = np.real(eigenvalues)
    eigenvectors = np.real(eigenvectors)

    print(f"For n = {n}, the values of lambda for which there are no solutions are:")

    lambda_values = []
    # 3. For each eigenvalue mu, find lambda = 1/mu and check the eigenvector condition
    for i in range(n):
        mu = eigenvalues[i]
        # The eigenvector from np.linalg.eig is a column vector
        c = eigenvectors[:, i]
        
        # Check if mu is non-zero to avoid division by zero
        if abs(mu) < 1e-9:
            # This case typically does not lead to a finite lambda
            continue
            
        # The condition for no solution is that the first component of the eigenvector is non-zero
        # We use a small tolerance for the check
        if abs(c[0]) > 1e-9:
            lambda_val = 1.0 / mu
            lambda_values.append(lambda_val)
            # The original equation format requested in the prompt
            print(f"λ = {lambda_val:.8f}")
            # Let's show the numbers in the original equation for clarity, setting n and λ
            # u(x) = 1 + {λ} ∫₀¹ dy (xⁿ - yⁿ) / (x - y) u(y)
            print(f"The equation with n={n} and this λ is:")
            print(f"u(x) = 1 + ({lambda_val:.8f}) * integral from 0 to 1 of (x^{n} - y^{n})/(x - y) * u(y) dy\n")

    if not lambda_values:
        print("No specific lambda values found under the given criteria.")


# --- Main execution ---
# The user can specify the value of n here. We will use n=3 as an example.
n_value = 3
solve_for_lambda(n_value)
