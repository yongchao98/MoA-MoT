import math

def solve_series():
    """
    Calculates the value of the series sum_{n=2 to inf} (-1)^n * S_n / n
    based on the derived formula: 1 - pi^2 / 12 + (ln(2))^2 / 2.
    """
    
    # Define the mathematical constants
    pi = math.pi
    ln2 = math.log(2)
    
    # Define the numbers in the final equation
    num1 = 1
    num2_numerator = pi**2
    num2_denominator = 12
    num3_numerator = ln2**2
    num3_denominator = 2
    
    # Calculate the terms
    term2 = num2_numerator / num2_denominator
    term3 = num3_numerator / num3_denominator
    
    # Calculate the final result
    result = num1 - term2 + term3
    
    # Print the equation with the values of its components
    print("The value of the sum is derived from the formula:")
    print(f"S = {num1} - (pi^2 / {num2_denominator}) + (ln(2)^2 / {num3_denominator})")
    print("\nSubstituting the values:")
    print(f"pi^2 = {num2_numerator}")
    print(f"ln(2)^2 = {num3_numerator}")
    print(f"S = {num1} - ({num2_numerator} / {num2_denominator}) + ({num3_numerator} / {num3_denominator})")
    print(f"S = {num1} - {term2} + {term3}")
    
    # Print the final result
    print("\nThe final value of the sum is:")
    print(result)

solve_series()