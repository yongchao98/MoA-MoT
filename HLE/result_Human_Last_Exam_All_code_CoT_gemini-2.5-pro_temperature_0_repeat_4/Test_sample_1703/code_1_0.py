import math

def solve_triangle_area():
    """
    This function calculates and prints the formula for the area of triangle T(t).
    
    The derivation is as follows:
    1.  The side length of the hexagon is R = 10.
    2.  The rotation of the hexagon does not affect the triangle's area.
    3.  The three lines containing the alternating sides of the hexagon form a large equilateral triangle
        with side length 3*R = 30. Its area is (sqrt(3)/4) * 30^2 = 225*sqrt(3).
    4.  The triangle T(t) is formed by cutting off three smaller triangles from the corners of this large triangle.
    5.  At time t, the vertices of T(t) have moved a distance of t from the midpoints of the hexagon's sides.
    6.  The area of each of the three identical corner triangles is (sqrt(3)/4) * (225 - t^2).
    7.  The area of T(t) is the area of the large triangle minus the areas of the three corner triangles:
        Area(t) = 225*sqrt(3) - 3 * (sqrt(3)/4) * (225 - t^2)
               = (sqrt(3)/4) * (900 - 675 + 3*t^2)
               = (sqrt(3)/4) * (225 + 3*t^2)
               = (3*sqrt(3)/4) * (75 + t^2)
    """
    
    # The numbers in the final equation Area(t) = (3 * sqrt(3) / 4) * (t^2 + 75)
    numerator_factor = 3
    sqrt_base = 3
    denominator = 4
    exponent = 2
    constant_term = 75
    
    # Print the final equation with each number explicitly shown.
    print("The area of the triangle T(t) as a function of time t is given by the formula:")
    print(f"Area(t) = ({numerator_factor} * sqrt({sqrt_base}) / {denominator}) * (t^{exponent} + {constant_term})")

solve_triangle_area()