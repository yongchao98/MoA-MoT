import math

def calculate_P(n):
    """
    Calculates the value of P(n) using the derived closed-form formula.
    
    The formula is P(n) = (2*pi)^(n*(n-1)/4) / sqrt(n!).
    
    Args:
        n (int): A positive integer.
        
    Returns:
        float: The calculated value of P(n).
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: Please provide a positive integer for n.")
        return

    # Calculate the components of the formula
    exponent = n * (n - 1) / 4
    try:
        factorial_val = math.factorial(n)
    except ValueError:
        print(f"Error: n={n} is too large to compute its factorial.")
        return

    # Print the formula with the specific numbers for the given n
    print(f"For n = {n}, the closed-form formula is:")
    print(f"P({n}) = (2 * pi)^({exponent}) / sqrt({factorial_val})")

    # Calculate the final result
    numerator = math.pow(2 * math.pi, exponent)
    denominator = math.sqrt(factorial_val)
    
    if denominator == 0:
        print("Error: Cannot divide by zero (denominator is zero).")
        return

    result = numerator / denominator
    
    # Print the final result
    print(f"The value is: P({n}) = {result}")
    return result

if __name__ == '__main__':
    # Example usage:
    try:
        n_input = int(input("Enter a positive integer n: "))
        calculate_P(n_input)
    except (ValueError, TypeError):
        print("Invalid input. Please enter a valid integer.")
