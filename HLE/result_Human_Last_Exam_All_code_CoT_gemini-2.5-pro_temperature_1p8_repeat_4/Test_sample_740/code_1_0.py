import numpy as np

def solve_lambda(n):
    """
    Calculates the values of lambda for which the integral equation has no solution.

    Args:
        n (int): The exponent in the integral equation kernel.
    """
    if not isinstance(n, int) or n < 1:
        print("n must be a positive integer.")
        return

    # Create the matrix A
    A = np.zeros((n, n))
    for j in range(n):
        for i in range(n):
            A[j, i] = 1.0 / (j - i + n)

    # Create the vector b
    b = np.zeros(n)
    for j in range(n):
        b[j] = 1.0 / (j + 1)
        
    print(f"For n = {n}, the system of linear equations is (I - λA)c = b")
    print("\nMatrix A:")
    # np.set_printoptions(precision=4, suppress=True) # Optional: for cleaner printing
    print(A)
    
    print("\nVector b:")
    print(b)
    
    # Calculate the eigenvalues of A
    try:
        eigenvalues = np.linalg.eigvals(A)
    except np.linalg.LinAlgError:
        print("\nCould not compute eigenvalues.")
        return
        
    # lambda values are the reciprocals of the eigenvalues
    lambda_values = 1.0 / eigenvalues

    print(f"\nThe values of λ for which the equation has no solution for n={n} are:")
    for val in lambda_values:
        # Checking if the imaginary part is negligible
        if np.iscomplex(val) and abs(val.imag) < 1e-9:
             print(f"{val.real:.8f}")
        else:
             print(f"{val:.8f}")


# Let's solve for n=3 as an example.
n_example = 3
solve_lambda(n_example)

# Let's demonstrate n=2 as well to compare with analytical results (-6 ± 4√3)
# -6 + 4*sqrt(3) ≈ 0.9282
# -6 - 4*sqrt(3) ≈ -12.9282
print("\n" + "="*30)
n_example = 2
solve_lambda(n_example)