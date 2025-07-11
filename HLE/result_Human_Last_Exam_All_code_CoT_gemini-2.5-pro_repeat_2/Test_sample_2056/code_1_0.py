import numpy as np
import sys

def calculate_l_k(n):
    """
    Calculates the value of l_k(n) based on the derived formula.

    The formula is:
    l_k(n) = (1/2)*ln(n+1) - 2*k^2 + k^2/n + (n-1)*ln(k)
    where k = ln(sqrt(2) + 1).

    Args:
        n (int): The dimension, must be >= 3.
    
    Returns:
        float: The calculated value of l_k(n).
    """
    # The constant k is defined as ln(sqrt(2) + 1)
    k = np.log(np.sqrt(2) + 1)

    # The derived formula for l_k(n) can be written as:
    # l_k(n) = c1 * log(n + c2) + c3 * k^2 + c4 * k^2 / n + (c5 * n + c6) * log(k)
    # where the numbers (coefficients) in the equation are:
    c1 = 0.5
    c2 = 1.0
    c3 = -2.0
    c4 = 1.0
    c5 = 1.0
    c6 = -1.0

    print("--- Calculating l_k(n) ---")
    print(f"The exact formula is: 0.5*ln(n+1) - 2*k^2 + k^2/n + (n-1)*ln(k)")
    print(f"Using n = {n} and k = {k:.8f}")

    # Calculate each term of the formula for the given n
    term1 = c1 * np.log(n + c2)
    term2 = c3 * k**2
    term3 = c4 * k**2 / n
    term4 = (c5 * n + c6) * np.log(k)
    
    result = term1 + term2 + term3 + term4

    print("\n--- Breakdown of Terms ---")
    print(f"Term 1 (0.5 * log(n+1)): {term1:.8f}")
    print(f"Term 2 (-2 * k^2): {term2:.8f}")
    print(f"Term 3 (k^2 / n): {term3:.8f}")
    print(f"Term 4 ((n-1) * log(k)): {term4:.8f}")
    
    print("\n--- Final Result ---")
    print(f"l_k({n}) = {result:.8f}")
    
    return result

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: python {sys.argv[0]} <n>")
        print("Description: Calculates l_k(n) for a given integer n >= 3.")
        sys.exit(1)
    try:
        n_val = int(sys.argv[1])
        if n_val < 3:
            print("Error: Input 'n' must be an integer greater than or equal to 3.")
            sys.exit(1)
        calculate_l_k(n_val)
    except ValueError:
        print("Error: Input 'n' must be a valid integer.")
        sys.exit(1)
