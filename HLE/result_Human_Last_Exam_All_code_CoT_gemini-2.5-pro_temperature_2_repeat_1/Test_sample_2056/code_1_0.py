import numpy as np

def calculate_l_k(n):
    """
    Calculates the exact value of l_k(n) for a given integer n >= 3.
    The function implements the derived analytical formula:
    l_k(n) = (n-1)*ln(k) + 0.5*ln(n+1) - k^2*(2 - 1/n)
    where k = ln(sqrt(2) + 1).

    Args:
        n (int): An integer, must be >= 3.
    """
    if not isinstance(n, int) or n < 3:
        print("Error: n must be an integer greater than or equal to 3.")
        return

    # Define the constant k as specified in the problem
    k = np.log(np.sqrt(2) + 1)

    # The derived analytical formula for l_k(n) consists of three main parts.
    # Term 1: (n-1) * ln(k)
    term1 = (n - 1) * np.log(k)

    # Term 2: (1/2) * ln(n+1)
    term2 = 0.5 * np.log(n + 1)

    # Term 3: -k^2 * (2 - 1/n)
    term3 = -k**2 * (2 - 1/n)

    # The final result is the sum of these three terms.
    result = term1 + term2 + term3
    
    # Print the breakdown of the calculation as requested
    print(f"Calculating l_k(n) for n = {n}:")
    print(f"The formula is: l_k(n) = (n-1)*ln(k) + 0.5*ln(n+1) - k^2*(2-1/n)")
    print(f"where k = ln(sqrt(2)+1) approx {k:.6f}\n")
    print("The numbers in the final equation are:")
    print(f"Term 1: ({n}-1)*ln({k:.6f}) = {term1:.6f}")
    print(f"Term 2: 0.5*ln({n}+1) = {term2:.6f}")
    print(f"Term 3: -({k:.6f})^2 * (2 - 1/{n}) = {term3:.6f}")
    print("-" * 30)
    print(f"l_k({n}) = {term1:.6f} + {term2:.6f} + {term3:.6f} = {result:.6f}")


# Example usage for n=5
n_value = 5
calculate_l_k(n_value)