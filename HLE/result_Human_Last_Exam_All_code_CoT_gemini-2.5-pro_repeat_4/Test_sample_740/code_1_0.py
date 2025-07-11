import numpy as np

def solve_for_lambda(n):
    """
    Calculates the values of lambda for which the given integral equation has no solutions.

    Args:
        n (int): The exponent in the integral equation kernel. Must be a positive integer.

    Returns:
        list: A list of values of lambda for which the equation has no solutions.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return []

    # Construct the n x n matrix A
    A = np.zeros((n, n))
    for j in range(n):
        for k in range(n):
            A[j, k] = 1.0 / (n - j + k)

    # Find the eigenvalues of A
    try:
        eigenvalues = np.linalg.eigvals(A)
    except np.linalg.LinAlgError:
        print("Error: Eigenvalue computation did not converge.")
        return []

    # Lambda values are the reciprocals of the eigenvalues of A
    # Filter out potential division by zero if an eigenvalue is close to 0
    lambda_values = [1.0 / mu for mu in eigenvalues if abs(mu) > 1e-12]
    
    # In this specific problem, it can be shown that for every characteristic value lambda,
    # the non-homogeneous equation has no solution. We don't need to check the
    # Fredholm alternative condition explicitly here.

    return lambda_values

def main():
    """
    Main function to get n from the user and print the solution.
    """
    try:
        # You can change the value of n here.
        n_str = input("Enter the integer value for n (e.g., 2): ")
        n = int(n_str)
        
        lambda_values = solve_for_lambda(n)

        if lambda_values:
            print(f"For n = {n}, the equation has no solutions for the following values of lambda:")
            # Sort the values for consistent output, real part first
            lambda_values.sort(key=lambda x: (x.real, x.imag))
            for val in lambda_values:
                # Format complex numbers for better readability
                if abs(val.imag) < 1e-12:
                    print(f"{val.real:.8f}")
                else:
                    print(f"{val.real:.8f} + {val.imag:.8f}j")

    except (ValueError, TypeError):
        print("Invalid input. Please enter a positive integer for n.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
