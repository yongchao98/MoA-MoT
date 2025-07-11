import math

def solve_for_p_n(n):
    """
    Calculates P(n) using the derived closed-form formula and prints the details.
    
    Args:
        n: A positive integer.
    """
    if not isinstance(n, int) or n <= 0:
        print("Input must be a positive integer.")
        return

    # The derived closed-form formula is P(n) = (2 * pi)^(n * (n - 1) / 4) / sqrt(n!)
    
    # Calculate components from the formula
    num_2 = 2
    num_pi_str = "pi" # For display
    num_n = n
    num_n_minus_1 = n - 1
    num_4 = 4
    
    # Use try-except for factorial to handle potential display issues with large n
    try:
        num_n_factorial = math.factorial(n)
        num_n_factorial_str = str(num_n_factorial)
    except (ValueError, OverflowError):
        num_n_factorial_str = f"{n}!" # Display as n! for very large n

    # To satisfy the instruction "output each number in the final equation!",
    # we print the formula with all numbers substituted for the given n.
    equation_str = f"P({num_n}) = ({num_2} * {num_pi_str})^({num_n} * {num_n_minus_1} / {num_4}) / sqrt({num_n_factorial_str})"
    print(f"The final equation for n={n} is:")
    print(equation_str)

    # Calculation using logarithms is more stable for larger numbers to avoid overflow/underflow
    exponent = n * (n - 1) / 4.0
    # log(P(n)) = exponent * log(2*pi) - 0.5 * log(n!)
    # log(n!) can be accurately calculated using math.lgamma(n + 1)
    log_p_n = exponent * math.log(2 * math.pi) - 0.5 * math.lgamma(n + 1)
    
    result = math.exp(log_p_n)
    
    print(f"\nThe calculated value is: P({n}) = {result}")


# Example calculation for n = 10
solve_for_p_n(10)