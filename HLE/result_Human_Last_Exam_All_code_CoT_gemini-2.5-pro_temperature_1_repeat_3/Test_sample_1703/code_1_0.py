import math

def print_area_formula():
    """
    This function prints the derived formula for the area of the triangle T(t).
    
    The derivation shows that the area A(t) is independent of the hexagon's rotation.
    Due to symmetry, the triangle remains equilateral.
    The area of an equilateral triangle is A = (sqrt(3)/4) * L^2, where L is the side length.
    
    By setting up a coordinate system for the hexagon with radius R=10 (and side length s=10),
    we can find the squared side length L(t)^2 of the triangle as a function of time t.
    The vertices of the triangle move with speed v=1 from the midpoints of alternating sides.
    
    The calculation yields the squared side length:
    L(t)^2 = 225 + 3 * t^2
    
    Substituting this into the area formula gives:
    A(t) = (sqrt(3) / 4) * (225 + 3 * t^2)
    
    This code will print the formula with its numerical constants.
    """
    
    # The constants in the final formula A(t) = (sqrt(c1) / c2) * (c3 + c4 * t^2)
    c1 = 3
    c2 = 4
    c3 = 225
    c4 = 3
    
    print("The area of the triangle T(t) as a function of time t is given by the formula:")
    # The following line prints the final equation with all the required numbers.
    print(f"A(t) = (sqrt({c1}) / {c2}) * ({c3} + {c4} * t^2)")

if __name__ == "__main__":
    print_area_formula()