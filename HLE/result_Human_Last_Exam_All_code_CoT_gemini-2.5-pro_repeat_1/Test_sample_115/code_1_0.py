import math

def solve_alpha():
    """
    Calculates the value of alpha based on the derivation.
    The growth rate of c_n is Theta(alpha^n), and we determined alpha to be the golden ratio phi.
    alpha = phi = (1 + sqrt(5)) / 2
    """
    
    # Define the numbers in the final equation
    num1 = 1
    num5 = 5
    num2 = 2
    
    # Calculate alpha
    alpha = (num1 + math.sqrt(num5)) / num2
    
    # Output the steps of the calculation and the final value
    print("The value of alpha is determined to be the golden ratio, phi.")
    print(f"The formula for alpha is: ({num1} + sqrt({num5})) / {num2}")
    print(f"The calculated value of alpha is: {alpha}")

solve_alpha()