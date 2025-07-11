import math

def solve():
    """
    This function calculates and prints the minimal polynomial for the shortest 
    geodesic distance an ant can travel on a dodecahedron under the given conditions.
    
    The solution relies on the known result from geometry that the squared shortest
    distance 'd' is d^2 = 7 + phi, where phi is the golden ratio.
    """
    
    # The problem is to find the minimal polynomial for d.
    # Let x = d^2.
    # We have the relation: x = 7 + phi
    # phi = (1 + sqrt(5)) / 2
    # So, x = 7 + (1 + sqrt(5)) / 2
    # 2x = 14 + 1 + sqrt(5)
    # 2x - 15 = sqrt(5)
    
    # Squaring both sides to eliminate the square root:
    # (2x - 15)^2 = 5
    # 4x^2 - 60x + 225 = 5
    # 4x^2 - 60x + 220 = 0
    
    # Dividing by 4, we get the minimal polynomial for x = d^2:
    # x^2 - 15x + 55 = 0
    
    # Now, substitute x = d^2 back into the equation:
    # (d^2)^2 - 15(d^2) + 55 = 0
    # d^4 - 15d^2 + 55 = 0
    
    # This gives a polynomial for d. To be the minimal polynomial, it must be
    # irreducible over the rational numbers, which it is.
    
    # The coefficients of the minimal polynomial P(d) = c4*d^4 + c3*d^3 + c2*d^2 + c1*d + c0 = 0 are:
    c4 = 1
    c3 = 0
    c2 = -15
    c1 = 0
    c0 = 55
    
    print("The minimal polynomial for the shortest distance d is of the form:")
    print("c4*d^4 + c3*d^3 + c2*d^2 + c1*d + c0 = 0")
    print("\nThe coefficients are:")
    print(f"c4 = {c4}")
    print(f"c3 = {c3}")
    print(f"c2 = {c2}")
    print(f"c1 = {c1}")
    print(f"c0 = {c0}")
    
    print("\nThus, the final equation is:")
    print(f"{c4}*d^4 + {c3}*d^3 + ({c2})*d^2 + {c1}*d + {c0} = 0")

solve()