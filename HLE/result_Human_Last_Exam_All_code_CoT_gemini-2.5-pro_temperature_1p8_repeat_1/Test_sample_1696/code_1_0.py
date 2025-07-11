import math

def calculate_P(n):
    """
    Calculates the value of P(n) using the derived closed-form formula.

    The closed-form formula for P(n) is:
    P(n) = (2 * pi)^(n * (n - 1) / 4) / sqrt(n!)
    
    Args:
        n (int): A positive integer.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: Please provide a positive integer for n.")
        return

    print(f"Calculating P(n) for n = {n}")
    print(f"The formula is: P(n) = (2 * pi)^(n * (n - 1) / 4) / sqrt(n!)")
    
    # Calculate the exponent for the numerator
    numerator_exponent = n * (n - 1) / 4
    
    # Calculate the numerator base
    base = 2 * math.pi
    
    # Calculate factorial of n for the denominator
    n_factorial = math.factorial(n)
    
    print(f"\nStep 1: Calculate the exponent of the numerator.")
    print(f"Exponent = {n} * ({n} - 1) / 4 = {numerator_exponent}")
    
    print(f"\nStep 2: Calculate the numerator: (2 * pi)^({numerator_exponent}).")
    numerator = math.pow(base, numerator_exponent)
    print(f"Numerator = ({base:.4f})^{numerator_exponent} = {numerator}")

    print(f"\nStep 3: Calculate the denominator: sqrt({n}!).")
    denominator = math.sqrt(n_factorial)
    print(f"Denominator = sqrt({n_factorial}) = {denominator:.4f}")
    
    # Final calculation
    if denominator == 0:
        result = "undefined (division by zero)"
    else:
        result = numerator / denominator

    print(f"\nFinal Result:")
    print(f"P({n}) = {numerator} / {denominator:.4f}")
    print(f"P({n}) = {result}")

if __name__ == '__main__':
    # You can change this value to calculate P(n) for any positive integer n.
    n_value = 10
    calculate_P(n_value)