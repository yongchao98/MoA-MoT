import math

def solve_for_p_n(n):
    """
    Calculates the value of P(n) using its derived closed-form formula.

    The closed form for P(n) is (2*pi)^(n*(n-1)/4) / sqrt(n!).
    This function computes the result for a given positive integer n and
    prints the components of the formula as well as the final value.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: Input 'n' must be a positive integer.")
        return

    # Calculate components of the formula
    n_minus_1 = n - 1
    exponent = n * n_minus_1 / 4
    
    try:
        n_factorial = math.factorial(n)
    except ValueError:
        print(f"Error: Cannot calculate factorial for n={n}.")
        return

    # Calculate the final result
    # (2*pi)^exponent
    val_numerator = (2 * math.pi) ** exponent
    # sqrt(n!)
    val_denominator = math.sqrt(n_factorial)
    
    if val_denominator == 0:
        print("Error: Division by zero (sqrt(n!) is zero).")
        return
        
    result = val_numerator / val_denominator
    
    # Output the formula with the numbers for the given n
    print("The derived closed-form formula for P(n) is: (2 * pi)^(n * (n - 1) / 4) / sqrt(n!)")
    print("-" * 30)
    print(f"For n = {n}:")
    print(f"The equation is: P({n}) = (2 * pi)^({n} * ({n} - 1) / 4) / sqrt({n}!)")
    print(f"                 P({n}) = (2 * {math.pi:.4f})^({exponent}) / sqrt({n_factorial})")
    print("-" * 30)
    print("Component values:")
    print(f"Value of n: {n}")
    print(f"Value of pi: {math.pi}")
    print(f"Value of n-1: {n_minus_1}")
    print(f"Value of the exponent (n * (n - 1) / 4): {exponent}")
    print(f"Value of n!: {n_factorial}")
    print("-" * 30)
    print(f"Final calculated value P({n}) = {result}")

# Set the value of n for the calculation.
# You can change this value to any positive integer.
n_input = 5
solve_for_p_n(n_input)