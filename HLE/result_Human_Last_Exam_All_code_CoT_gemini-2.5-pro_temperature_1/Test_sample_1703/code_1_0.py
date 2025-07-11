import math

def solve_triangle_area():
    """
    This function calculates and prints the formula for the area of triangle T(t).
    """
    
    # Define the constants from the problem statement.
    R = 10  # Radius of the circle
    v = 1   # Speed of the vertices
    
    # The side length of the regular hexagon is equal to the radius.
    s = R
    
    # As derived in the plan, the area of the triangle T(t) is determined by its
    # squared side length, a(t)^2. The rotation of the hexagon does not affect the area.
    
    # The initial squared side length of the equilateral triangle at t=0, when vertices
    # are at the midpoints of alternating sides, is (1.5 * R)^2.
    initial_squared_side_length = (1.5 * R)**2
    
    # The change in the squared side length is due to the vertices moving with speed v=1.
    # The derived formula for the squared side length at time t is:
    # a(t)^2 = initial_squared_side_length + 3 * t^2
    # a(t)^2 = 225 + 3*t^2
    
    # The area of an equilateral triangle with side 'a' is (sqrt(3)/4) * a^2.
    # So, Area(t) = (sqrt(3)/4) * (225 + 3*t^2)
    # This can be expanded to: Area(t) = (225 * sqrt(3) / 4) + (3 * sqrt(3) / 4) * t^2
    
    # The problem asks to output each number in the final equation.
    # We will construct and print the string representing the formula.
    
    # Coefficients of the formula Area(t) = C1 + C2 * t^2
    c1_num = int(initial_squared_side_length)
    c1_den = 4
    
    c2_num = 3
    c2_den = 4
    
    # Print the final equation for the area as a function of t.
    print("The area of the triangle T(t) as a function of time is given by the equation:")
    print(f"Area(t) = ({c1_num} * sqrt(3) / {c1_den}) + ({c2_num} * sqrt(3) / {c2_den}) * t^2")

    # To show the numerical value of the coefficients:
    coeff1 = (c1_num * math.sqrt(3)) / c1_den
    coeff2 = (c2_num * math.sqrt(3)) / c2_den
    print("\nNumerically, the function is approximately:")
    print(f"Area(t) = {coeff1:.4f} + {coeff2:.4f} * t^2")

solve_triangle_area()