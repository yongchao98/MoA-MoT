import math

def solve_p_n(n):
    """
    Calculates the value of P(n) using the derived closed-form formula.
    P(n) = (2*pi)^(n*(n-1)/4) / sqrt(n!)
    """
    if not isinstance(n, int) or n <= 0:
        print("Please provide a positive integer for n.")
        return

    # Calculate components of the formula
    n_factorial = math.factorial(n)
    sqrt_n_factorial = math.sqrt(n_factorial)
    
    exponent = n * (n - 1) / 4.0
    two_pi = 2 * math.pi
    numerator = math.pow(two_pi, exponent)
    
    result = numerator / sqrt_n_factorial

    # Print the formula with numbers
    print(f"The closed-form formula for P(n) is (2*pi)**(n*(n-1)/4) / sqrt(n!)")
    print(f"For n = {n}:")
    print(f"P({n}) = (2 * {math.pi})**({n} * ({n}-1) / 4) / sqrt({n}!)")
    print(f"P({n}) = ({two_pi})**({exponent}) / sqrt({n_factorial})")
    print(f"P({n}) = {numerator} / {sqrt_n_factorial}")
    print(f"P({n}) = {result}")

# Example usage with n=10 as a demonstration
# The user can change this value.
n = 10
solve_p_n(n)