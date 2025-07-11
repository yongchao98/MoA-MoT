import math

def calculate_limit():
    """
    Calculates the limit of n*P(n) as n -> infinity.
    The derived analytical limit is 2*sqrt(3)/pi.
    """
    numerator_val = 2
    sqrt_val = 3
    
    # Calculate the limit value
    result = numerator_val * math.sqrt(sqrt_val) / math.pi
    
    # Print the equation and the result
    print(f"The limit is given by the expression: ({numerator_val} * sqrt({sqrt_val})) / pi")
    print(f"The numerical value of the limit is: {result}")

if __name__ == "__main__":
    calculate_limit()
