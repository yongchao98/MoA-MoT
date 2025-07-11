import math

def calculate_deflection():
    """
    Calculates the membrane's deflection at x = 0.
    The deflection y(0) is found by solving the Lagrange differential equation,
    which leads to the final expression: y(0) = ((81/16)^(4/5)) / 3.
    """
    
    # Define the numbers in the final equation
    numerator_base_num = 81
    numerator_base_den = 16
    exponent_num = 4
    exponent_den = 5
    denominator = 3
    
    # Print the final equation
    print(f"The final equation for the deflection y(0) is: y(0) = (({numerator_base_num} / {numerator_base_den})^({exponent_num}/{exponent_den})) / {denominator}")
    
    # Calculate the value
    base = numerator_base_num / numerator_base_den
    exponent = exponent_num / exponent_den
    
    y_at_0 = (base ** exponent) / denominator
    
    # Print the calculated result
    print("\nThe calculated value for the deflection at x = 0 is:")
    print(y_at_0)

calculate_deflection()