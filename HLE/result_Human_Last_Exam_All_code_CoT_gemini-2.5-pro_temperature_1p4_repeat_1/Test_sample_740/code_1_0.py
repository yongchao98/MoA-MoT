import numpy as np

def find_lambda_no_solution(n):
    """
    For a given integer n, this function calculates the values of λ for which the
    integral equation u(x) = 1 + λ ∫[0,1] dy (xⁿ - yⁿ)/(x-y) u(y) has no solutions.

    The method involves converting the integral equation into a system of linear
    equations (I - λA)c = b. The equation has no solutions for the characteristic
    values of λ, which are the reciprocals of the eigenvalues of the matrix A.

    Args:
        n (int): The positive integer exponent in the integral equation.
    
    Returns:
        np.ndarray: An array containing the n values of λ.
    """
    if not isinstance(n, int) or n <= 0:
        raise ValueError("n must be a positive integer.")

    # Construct the n x n matrix A, where A_jk = 1 / (j + n - k)
    # using 0-based indexing for row j and column k.
    A = np.fromfunction(lambda j, k: 1.0 / (j - k + n), (n, n), dtype=int)

    # The values of 1/λ are the eigenvalues of matrix A.
    eigenvalues = np.linalg.eigvals(A)

    # The eigenvalues of A for this problem are real. We can sort them.
    # We take the real part to discard negligible imaginary parts from numerical error.
    sorted_eigenvalues = np.sort(eigenvalues.real)
    
    # The values of λ for which there are no solutions are the reciprocals
    # of the eigenvalues.
    lambda_values = 1.0 / sorted_eigenvalues
    
    return lambda_values

def main():
    """
    Main function to demonstrate the solution for a sample value of n.
    """
    try:
        # The problem statement has n as a parameter. We will demonstrate for n=3.
        n_example = 3
        print(f"For n = {n_example}, the values of λ for which the equation has no solution are:")
        
        lambda_values = find_lambda_no_solution(n_example)
        
        # The problem asks to output each number in the final equation.
        # We will print each calculated λ value.
        for i, val in enumerate(lambda_values):
            print(f"λ_{i+1} = {val:.8f}")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()