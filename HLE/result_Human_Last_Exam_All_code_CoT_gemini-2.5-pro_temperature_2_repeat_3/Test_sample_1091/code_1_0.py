import math

def calculate_limit():
    """
    Calculates the limit of n*P(n) as n -> infinity.
    The final expression for the limit is (2 * sqrt(3)) / pi.
    This script calculates the numerical value and prints the components of the expression.
    """
    
    # Define the components of the final expression
    val_2 = 2
    val_sqrt3 = math.sqrt(3)
    val_pi = math.pi
    
    # Calculate the final result
    result = (val_2 * val_sqrt3) / val_pi
    
    # Print the equation with substituted values and the final result
    print("The analytical expression for the limit is (2 * sqrt(3)) / pi")
    print("Calculating the numerical value:")
    print(f"({val_2} * {val_sqrt3:.5f}) / {val_pi:.5f} = {result:.5f}")

calculate_limit()