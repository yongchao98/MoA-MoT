import numpy as np

def find_lambda_for_no_solution():
    """
    Finds the values of λ for which the given integral equation has no solution.
    The integer n is obtained from user input.
    """
    try:
        n_input = input("Enter the integer value for n (e.g., 2): ")
        n = int(n_input)
        if n <= 0:
            print("Error: n must be a positive integer.")
            return
    except ValueError:
        print("Error: Invalid input. Please enter an integer.")
        return

    # Step 1: Construct the matrix A and vector b
    # A_ij = 1/(i-j+n) for i,j in {0, ..., n-1}
    # b_i = 1/(i+1) for i in {0, ..., n-1}
    A = np.zeros((n, n), dtype=float)
    b = np.zeros(n, dtype=float)
    for i in range(n):
        b[i] = 1.0 / (i + 1)
        for j in range(n):
            A[i, j] = 1.0 / (i - j + n)

    # Step 2: Compute eigenvalues (mu) and eigenvectors (v) of A transpose
    try:
        # A.T @ v = mu * v
        eigenvalues_mu, eigenvectors_v = np.linalg.eig(A.T)
    except np.linalg.LinAlgError:
        print(f"Could not compute eigenvalues for the matrix with n={n}.")
        return

    # Step 3: Check the condition for no solution for each eigenvalue
    no_solution_lambdas = []
    for i in range(len(eigenvalues_mu)):
        mu = eigenvalues_mu[i]
        v = eigenvectors_v[:, i]

        # We must avoid division by zero if an eigenvalue is zero.
        if np.isclose(mu, 0):
            continue

        # Check if b is orthogonal to the eigenvector v
        # v.conj().T @ b is the dot product for complex vectors.
        dot_product = np.dot(v.conj(), b)

        # If the dot product is not zero, then no solution exists.
        if not np.isclose(dot_product, 0):
            lambda_val = 1.0 / mu
            no_solution_lambdas.append(lambda_val)

    # Step 4: Print the results
    print("-" * 50)
    print(f"For n = {n}, the values of λ for which the equation has no solutions are:")
    if not no_solution_lambdas:
        print("None found. A solution exists for all characteristic values of λ.")
    else:
        # Sort for consistent output. Real and imag parts are sorted separately for complex numbers.
        no_solution_lambdas.sort(key=lambda x: (x.real, x.imag))
        for lam in no_solution_lambdas:
            print(f"λ = {lam}")
            
    # For n=2, we can demonstrate the symbolic equation.
    if n == 2:
        print("\nFor the specific case n=2, these values are the roots of the quadratic equation:")
        a, b, c = 1, 12, -12
        print(f"{a}λ² + {b}λ + {c} = 0")
        print("The symbolic solution is λ = -6 ± 4√3.")

if __name__ == '__main__':
    find_lambda_for_no_solution()
