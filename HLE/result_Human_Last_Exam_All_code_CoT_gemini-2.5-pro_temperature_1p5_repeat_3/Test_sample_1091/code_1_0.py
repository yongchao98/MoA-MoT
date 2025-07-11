import math

def calculate_limit():
    """
    This function calculates the limit of n*P(n) as n -> infinity.
    The theoretical derivation shows that this limit is equal to (2 * sqrt(3)) / pi.
    """
    
    # The numbers in the final expression for the limit
    numerator_factor1 = 2
    numerator_factor2 = math.sqrt(3)
    denominator = math.pi
    
    # The final result
    result = (numerator_factor1 * numerator_factor2) / denominator
    
    # Printing the components of the equation and the final result
    print("The limit is given by the expression: (2 * sqrt(3)) / pi")
    print(f"Value of 2: {numerator_factor1}")
    print(f"Value of sqrt(3): {numerator_factor2}")
    print(f"Value of pi: {denominator}")
    print(f"The numerical value of the limit is: {result}")

calculate_limit()