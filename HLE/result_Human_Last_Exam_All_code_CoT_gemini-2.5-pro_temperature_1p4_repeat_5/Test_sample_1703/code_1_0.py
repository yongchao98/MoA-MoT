import math

def solve_triangle_area():
    """
    Calculates and prints the formula for the area of the triangle T(t).
    """
    # Step 1: Define constants from the problem description.
    # Radius of the circle C.
    R = 10
    # Speed of the triangle's vertices along the hexagon's sides.
    v = 1

    # Step 2: Determine the properties of the hexagon and the initial triangle T(0).
    # The side length of a regular hexagon inscribed in a circle is equal to the radius.
    hexagon_side_length = R
    
    # The vertices of T(0) are the midpoints of three alternating sides.
    # By geometric construction, this forms an equilateral triangle.
    # The side length of this initial triangle, L0, is sqrt(3) times the distance from the center
    # to a side, which is (sqrt(3)/2)*s. So, L0 = sqrt(3) * (sqrt(3)/2 * s) * sqrt(3) is wrong.
    # A direct calculation shows the side length L0 is 1.5 times the hexagon side length.
    L0 = 1.5 * hexagon_side_length
    
    # Step 3: Determine the formula for the squared side length of T(t).
    # Based on vector analysis of the vertices' motion, the squared side length L(t)^2 is:
    # L(t)^2 = L0^2 + 3 * (v*t)^2
    L0_sq = L0**2
    t_squared_coefficient = 3 * (v**2)

    # Step 4: Determine the formula for the area of T(t).
    # The area of an equilateral triangle is (sqrt(3)/4) * L^2.
    # Area(t) = (sqrt(3)/4) * (L0^2 + 3*(v*t)^2)
    # We will print this formula showing the calculated numbers.

    print("The area of the triangle T(t) as a function of time t is given by the formula:")
    print("Area(t) = (sqrt(3) / 4) * (L_squared)")
    print("where L_squared is the side length of the triangle squared.")
    print("\nBased on the problem's parameters:")
    print(f"L_squared(t) = {L0_sq} + {t_squared_coefficient}*t^2")
    
    print("\nSubstituting this into the area formula, we get the final equation:")
    # The instruction is to output each number in the final equation.
    # The final equation is Area(t) = (sqrt(3)/4) * (225 + 3*t^2)
    # The numbers are sqrt(3), 4, 225, 3.
    print(f"Area(t) = (sqrt(3) / 4) * ({int(L0_sq)} + {int(t_squared_coefficient)} * t^2)")
    
    # Expanding the formula
    constant_term_numerator = int(L0_sq)
    t_squared_term_numerator = int(t_squared_coefficient)
    denominator = 4
    
    print("\nExpanded form:")
    print(f"Area(t) = ({constant_term_numerator}*sqrt(3) / {denominator}) + ({t_squared_term_numerator}*sqrt(3) / {denominator}) * t^2")

solve_triangle_area()