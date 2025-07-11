import math

def solve_probability():
    """
    Calculates the approximate probability based on the derived formula.
    """
    # The target point is (3600, 0), so the radius R is 3600.
    R = 3600
    
    # The derived formula for the probability is P = 2 / (pi * R)
    numerator = 2
    denominator_R = 3600
    
    # Calculate the probability
    probability = numerator / (math.pi * denominator_R)
    
    # Print the equation as requested
    print(f"The probability P is calculated using the formula:")
    print(f"P = {numerator} / (pi * {denominator_R})")
    
    # Print the result
    print(f"\nCalculated probability: {probability}")
    
    # Print the approximate answer with two significant digits
    print(f"Approximate answer with two significant digits: {probability:.2g}")

solve_probability()