import math

def solve_p_n(n):
    """
    Calculates P(n) using the derived closed-form formula.
    
    The formula is P(n) = (2*pi)^(n*(n-1)/4) / sqrt(n!)
    
    Args:
        n (int): A positive integer.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: Please provide a positive integer for n.")
        return

    # Calculate the components of the formula
    exponent_numerator = n * (n - 1)
    exponent = exponent_numerator / 4
    n_factorial = math.factorial(n)
    sqrt_n_factorial = math.sqrt(n_factorial)
    
    # Calculate the final value of P(n)
    # Using logarithms for numerical stability with large n
    log_p_n = exponent * math.log(2 * math.pi) - 0.5 * math.log(n_factorial)
    p_n_value = math.exp(log_p_n)

    # Output the formula with the calculated numbers
    print(f"For n = {n}, the closed formula for P(n) is (2*pi)^(n*(n-1)/4) / sqrt(n!)")
    print("\nCalculating the components of the formula:")
    print(f"n = {n}")
    print(f"n-1 = {n-1}")
    print(f"The exponent is n*(n-1)/4 = {n}*({n-1})/4 = {exponent_numerator}/4 = {exponent}")
    print(f"n! = {n_factorial}")
    
    print(f"\nFinal Equation:")
    print(f"P({n}) = (2*pi)^({exponent}) / sqrt({n_factorial})")
    
    print(f"\nNumerical Value:")
    print(f"P({n}) = {p_n_value}")

if __name__ == '__main__':
    # Example usage with n=4
    n_value = 4
    solve_p_n(n_value)
    
    print("\n" + "="*20 + "\n")
    
    # Example usage with n=10
    n_value = 10
    solve_p_n(n_value)