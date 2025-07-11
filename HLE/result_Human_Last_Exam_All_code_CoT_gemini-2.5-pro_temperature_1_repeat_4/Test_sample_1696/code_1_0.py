import math

def calculate_P(n):
    """
    Calculates P(n) using the derived closed-form formula.
    P(n) = (2*pi)^(n*(n-1)/4) / sqrt(n!)
    
    Args:
        n (int): A positive integer.
        
    Returns:
        A tuple containing the result, the exponent, and n!
    """
    if not isinstance(n, int) or n <= 0:
        raise ValueError("n must be a positive integer")

    # Calculate components of the formula
    exponent = n * (n - 1) / 4
    try:
        n_factorial = math.factorial(n)
    except ValueError:
        # Handle cases where n is too large for standard float
        # This part of code will not be hit with reasonable n
        print("n is too large to compute factorial directly.")
        return None, None, None
    
    # Calculate P(n) using logarithms for numerical stability with large n
    log_p_n = exponent * math.log(2 * math.pi) - 0.5 * math.log(n_factorial)
    result = math.exp(log_p_n)
    
    return result, exponent, n_factorial

def main():
    """
    Main function to demonstrate the formula for a given n.
    """
    # We can choose any positive integer n for demonstration.
    # Let's use n=4 as it's a non-trivial case involving both product rules.
    n = 4
    
    try:
        p_n, exponent, n_factorial = calculate_P(n)

        # Print the final formula in a readable way
        print(f"The closed-form formula for P(n) is:")
        print(f"P(n) = (2 * pi)^(n * (n-1) / 4) / sqrt(n!)")
        
        print("\nSubstituting the value of n:")
        # As requested: "output each number in the final equation!"
        print(f"P({n}) = (2 * pi)^({n} * ({n}-1) / 4) / sqrt({n}!)")
        print(f"P({n}) = (2 * pi)^({int(exponent)}) / sqrt({n_factorial})")
        
        # Print the final numerical result
        print(f"\nThe calculated value of P({n}) is: {p_n}")

    except (ValueError, TypeError) as e:
        print(e)

if __name__ == "__main__":
    main()