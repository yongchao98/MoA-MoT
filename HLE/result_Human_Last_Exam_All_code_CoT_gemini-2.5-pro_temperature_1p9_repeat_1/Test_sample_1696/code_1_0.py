import math

def calculate_P(n):
    """
    Calculates the value of P(n) using the derived closed-form formula.
    
    The formula is P(n) = (2*pi)^(n*(n-1)/4) / sqrt(n!).
    To avoid overflow with large numbers, calculations are done using logarithms.
    log(P(n)) = (n*(n-1)/4) * log(2*pi) - 0.5 * log(n!)
    
    Args:
        n (int): A positive integer.
        
    Returns:
        float: The calculated value of P(n), or None if input is invalid.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: Please provide a positive integer for n.")
        return None

    print(f"Finding the value of P(n) for n = {n}")
    print("The closed-form formula is: P(n) = (2 * pi)^(n * (n - 1) / 4) / sqrt(n!)")
    print("-" * 30)

    # Numerator part
    exponent = n * (n - 1) / 4
    log_2pi = math.log(2 * math.pi)
    log_numerator = exponent * log_2pi
    
    # Denominator part (using lgamma for log(n!))
    # math.lgamma(n + 1) is log(n!)
    log_factorial_n = math.lgamma(n + 1)
    log_denominator = 0.5 * log_factorial_n
    
    # Final result in log scale, then converted back
    log_P_n = log_numerator - log_denominator
    result = math.exp(log_P_n)
    
    # Print the equation with substituted values
    print(f"Substituting n = {n}:")
    n_fact = math.factorial(n)
    
    print(f"P({n}) = (2 * {math.pi:.6f})^({n} * {n-1} / 4) / sqrt({n}!)")
    print(f"P({n}) = ({2 * math.pi:.6f})^{exponent} / sqrt({n_fact})")
    
    num_val = math.exp(log_numerator)
    den_val = math.exp(log_denominator)

    # For very large n, num_val or den_val could be inf, handle this for printing
    if math.isinf(num_val) or math.isinf(den_val):
        print("Intermediate values are too large to display directly.")
        print(f"log(Numerator) = {log_numerator:.6f}")
        print(f"log(Denominator) = {log_denominator:.6f}")
        print(f"log(P({n})) = {log_P_n:.6f}")
    else:
        print(f"P({n}) = {num_val:.6f} / {den_val:.6f}")
        
    print("-" * 30)
    print(f"The final value of P({n}) is: {result}")
    
    return result

if __name__ == '__main__':
    # Get user input for n
    try:
        n_input_str = input("Enter a positive integer n: ")
        n_input = int(n_input_str)
        calculate_P(n_input)
    except (ValueError, TypeError):
        print("Invalid input. Please enter a positive integer.")
