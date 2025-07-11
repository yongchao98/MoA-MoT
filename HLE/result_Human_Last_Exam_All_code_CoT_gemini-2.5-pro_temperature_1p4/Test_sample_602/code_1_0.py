import numpy as np

def calculate_l_n(n):
    """
    Calculates the value of the function l(n) for a given integer n >= 5.

    The function implements the derived exact formula:
    l(n) = (2/n^2) * [n^2 + 1 - (2n-1) * sqrt(n^2 - n + 1)]
    
    It prints the components of the formula and the final numerical result.
    """
    # Validate the input n
    if not isinstance(n, int) or n < 5:
        print("Error: n must be an integer greater than or equal to 5.")
        return

    # Calculate the components of the formula
    n_squared = n**2
    n_squared_plus_1 = n_squared + 1
    two_n_minus_1 = 2 * n - 1
    radicand = n_squared - n + 1
    sqrt_val = np.sqrt(radicand)

    # Calculate the final value using the formula
    term_in_brackets = n_squared_plus_1 - two_n_minus_1 * sqrt_val
    l_n_value = (2 / n_squared) * term_in_brackets

    # Output the results in a structured way
    print(f"Calculating l(n) for n = {n}")
    print(f"The formula is: l(n) = (2/n^2) * [n^2 + 1 - (2n - 1) * sqrt(n^2 - n + 1)]")
    print("\n--- Intermediate values ---")
    print(f"n^2              = {n_squared}")
    print(f"n^2 + 1          = {n_squared_plus_1}")
    print(f"2n - 1           = {two_n_minus_1}")
    print(f"n^2 - n + 1      = {radicand}")
    print(f"sqrt(n^2 - n + 1) = {sqrt_val:.8f}")
    
    print("\n--- Final equation using calculated values ---")
    print(f"l({n}) = (2/{n_squared}) * [{n_squared_plus_1} - {two_n_minus_1} * {sqrt_val:.8f}]")
    print(f"l({n}) = {l_n_value:.8f}")

# Example calculation for n=5, the smallest value in the domain.
# You can change this value to calculate l(n) for other integers >= 5.
n_val = 5
calculate_l_n(n_val)
