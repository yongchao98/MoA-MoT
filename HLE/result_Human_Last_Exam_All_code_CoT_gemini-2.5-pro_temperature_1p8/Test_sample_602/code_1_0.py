import math

def calculate_l(n: int):
    """
    Calculates the value of l(n) for n >= 5 based on the derived formula.
    
    The derived formula for l(n) is:
    l(n) = (2/n^2) * (n^2 + 1 - (2n-1)*sqrt(n^2 - n + 1))
    
    Args:
        n: An integer greater than or equal to 5.
        
    Returns:
        The float value of l(n).
    """
    if not isinstance(n, int) or n < 5:
        raise ValueError("Input must be an integer n >= 5.")

    n_squared = n * n
    
    # Components of the equation
    term_a = 2.0
    term_b = n_squared
    term_c = n_squared + 1.0
    term_d = 2.0 * n - 1.0
    sqrt_arg = n_squared - n + 1.0
    term_e = math.sqrt(sqrt_arg)
    
    # Printing the numbers in the final equation for the given n
    print(f"Calculating l(n) for n = {n}:")
    print(f"l({n}) = ({term_a} / {n}^2) * ({n}^2 + 1 - (2*{n}-1) * sqrt({n}^2 - {n} + 1))")
    print(f"l({n}) = ({term_a} / {term_b}) * ({term_c} - {term_d} * sqrt({sqrt_arg}))")
    print(f"l({n}) = ({term_a} / {term_b}) * ({term_c} - {term_d} * {term_e})")

    # Final calculation
    result = (term_a / term_b) * (term_c - term_d * term_e)
    
    print("\nFinal Result:")
    print(f"l({n}) = {result}")
    
    return result

if __name__ == '__main__':
    # You can change the value of n here.
    # The problem is defined for n >= 5.
    n_value = 10
    calculate_l(n_value)
