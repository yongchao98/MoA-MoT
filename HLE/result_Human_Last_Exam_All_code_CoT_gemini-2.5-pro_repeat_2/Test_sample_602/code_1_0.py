import numpy as np

def solve_for_n(n):
    """
    Solves the problem for a given integer n >= 5.
    """
    if not isinstance(n, int) or n < 5:
        raise ValueError("n must be an integer greater than or equal to 5.")

    # Step 1: Define matrices M and P
    a = np.sqrt(1 - (n - 1) / n**2)
    b = 1 / n
    M = np.full((n, n), b)
    np.fill_diagonal(M, a)

    # Step 2: Calculate f^{(3)}(P).
    # This is the inverse of P. From analytical derivation, we know its structure.
    X_prime = np.zeros((n, n))
    np.fill_diagonal(X_prime, 2)
    np.fill_diagonal(X_prime[1:], 1)
    np.fill_diagonal(X_prime[:, 1:], 1)

    # Step 3: Calculate f^{(2)}(X_prime), the projection onto the tangent space.
    # The projected matrix is Z = X_prime - M * mdiag(M^T * X_prime)
    lambda_diag_vector = np.diag(M.T @ X_prime)
    Lambda_matrix = np.diag(lambda_diag_vector)
    Z = X_prime - M @ Lambda_matrix

    # Step 4: Calculate f^{(1)}(Z), the sum of the first and last rows.
    result = np.sum(Z[0, :]) + np.sum(Z[n-1, :])
    
    return result

def main():
    """
    Main function to execute the solution.
    """
    # As per the prompt, n >= 5. Let's use n=5 as an example.
    n = 5
    l_n_value = solve_for_n(n)
    
    print(f"The numerical value of l({n}) is: {l_n_value}")

    # The prompt also asks to "output each number in the final equation!".
    # The derived analytical formula is:
    # l(n) = (2*n**2 + 2 - 2*(2*n - 1)*sqrt(n**2 - n + 1)) / n**2
    # We can rewrite this as: n^2 * l(n) = 2*n^2 + 2 - 2*(2*n - 1)*sqrt(1*n^2 - 1*n + 1)
    # The numbers in order of appearance are printed below.
    print("\nThe numbers from the final derived equation are:")
    print(2) # from 2*n^2
    print(2) # from n^2
    print(2) # from +2
    print(-2) # from -2*(...)
    print(2) # from 2*n
    print(-1) # from -1
    print(1) # from 1*n^2
    print(-1) # from -1*n
    print(1) # from +1

if __name__ == "__main__":
    main()
