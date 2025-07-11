import math

def solve_l(n: int):
    """
    Calculates the value of l(n) for n >= 5.

    The function follows the derived analytical formula:
    l(n) = (2/n^2) * (n^2 + 1 - (2n - 1) * sqrt(n^2 - n + 1))
    """
    if not isinstance(n, int) or n < 5:
        print("Error: Input must be an integer n >= 5.")
        return

    # Components of the derived formula
    n_squared = n * n
    term1 = n_squared + 1
    term2_coeff = 2 * n - 1
    sqrt_term_arg = n_squared - n + 1
    
    # Calculate the argument of the square root
    try:
        sqrt_term_val = math.sqrt(sqrt_term_arg)
    except ValueError:
        print(f"Error: Cannot compute square root of a negative number for n={n}.")
        return

    # Construct the equation string
    # This fulfills the requirement to "output each number in the final equation"
    equation_str = f"l({n}) = (2 / {n_squared}) * ({term1} - {term2_coeff} * sqrt({sqrt_term_arg}))"
    
    # Calculate the final value
    result = (2 / n_squared) * (term1 - term2_coeff * sqrt_term_val)

    print(f"For n = {n}:")
    print(f"The equation is: {equation_str}")
    print(f"The exact value is: l({n}) = {result}")

# Example usage with n=5 as required by the problem's domain N >= 5
# You can change this value to any integer >= 5.
n_value = 5
solve_l(n_value)

# Example usage with another value, e.g., n=10
n_value = 10
solve_l(n_value)