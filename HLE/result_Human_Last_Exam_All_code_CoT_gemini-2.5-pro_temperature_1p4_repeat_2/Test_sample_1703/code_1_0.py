import math

def solve_triangle_area():
    """
    This function calculates and prints the formula for the area of the triangle T(t).

    The derivation steps are:
    1. The rotation of the hexagon does not affect the area of the inscribed triangle.
    2. The triangle T(t) starts as an equilateral triangle with side length a_0 = 1.5 * S,
       where S=10 is the side length of the hexagon. So, a_0 = 15.
    3. The initial area A(0) is (sqrt(3)/4) * 15^2 = 225 * sqrt(3) / 4.
    4. The vertices move with speed v=1, so their displacement from the midpoint at time t is t.
    5. Due to symmetry, the triangle T(t) remains equilateral. Its side length squared a(t)^2 can be
       shown to be a(t)^2 = a_0^2 + 3*t^2 = 225 + 3*t^2.
    6. The area A(t) is (sqrt(3)/4) * a(t)^2.
    """
    
    # Constants from the derivation
    a0_squared = 225
    t_squared_coefficient = 3
    
    # The area formula is A(t) = (sqrt(3)/4) * (a0^2 + 3*t^2)
    # We will print this formula piece by piece as requested.
    
    print("The area of the triangle T(t) as a function of time t is given by the formula:")
    print("Area(t) = (sqrt(3) / 4) * (a_0^2 + k * t^2)")
    print("\nWhere the numbers in the equation are:")
    print(f"The constant for the square root of 3 is: {math.sqrt(3):.4f}...")
    print("The constant divisor is: 4")
    print(f"The initial side length squared, a_0^2, is: {a0_squared}")
    print(f"The coefficient for t^2, k, is: {t_squared_coefficient}")
    
    print("\nSo the final formula is:")
    print(f"Area(t) = (sqrt(3) / 4) * ({a0_squared} + {t_squared_coefficient} * t^2)")
    
    # This formula is valid for t in [0, 5]. For t > 5, the vertices would
    # reach the end of their respective sides and presumably reverse direction.
    # In that case, t in the formula would be replaced by a periodic function.
    # Assuming the question implies the formula for the initial phase of motion.

solve_triangle_area()