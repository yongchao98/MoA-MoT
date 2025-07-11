import math

def calculate_lower_bound():
    """
    Calculates the constant lower bound for d(t,x).
    
    The lower bound is the minimum value of the function m_-(u) for u in [0, 1].
    The function m_-(u) is the lower root of a quadratic derived from the PDE.
    m_-(u) = (3*u - 5*u^2 - u * sqrt(17*u^2 - 22*u + 9)) / 4
    
    The minimum of this function occurs at u=1. This script calculates m_-(1).
    """
    
    u = 1.0
    
    # Components of the equation for m_-(1)
    term1 = 3 * u
    term2 = -5 * u**2
    
    # Inside the square root
    sqrt_term_inside_1 = 17 * u**2
    sqrt_term_inside_2 = -22 * u
    sqrt_term_inside_3 = 9
    
    sqrt_val = math.sqrt(sqrt_term_inside_1 + sqrt_term_inside_2 + sqrt_term_inside_3)
    
    term3 = -u * sqrt_val
    
    numerator = term1 + term2 + term3
    denominator = 4.0
    
    lower_bound = numerator / denominator
    
    print("The lower bound is found by calculating the minimum of a function m_-(u) over u in [0,1].")
    print("The minimum occurs at u = 1.")
    print("\nThe equation is: m_-(1) = (3*1 - 5*1^2 - 1*sqrt(17*1^2 - 22*1 + 9)) / 4")
    print("\nCalculation steps:")
    print(f"Numerator part 1 (3*u): {term1}")
    print(f"Numerator part 2 (-5*u^2): {term2}")
    print(f"Value inside square root (17*u^2 - 22*u + 9): {sqrt_term_inside_1 + sqrt_term_inside_2 + sqrt_term_inside_3}")
    print(f"Square root value: {sqrt_val}")
    print(f"Numerator part 3 (-u*sqrt(...)): {term3}")
    print(f"Total numerator: {numerator}")
    print(f"Denominator: {denominator}")
    
    print("\nFinal constant lower bound:")
    print(lower_bound)

calculate_lower_bound()