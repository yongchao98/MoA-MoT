import math

def solve_triangle_area():
    """
    This function derives and prints the formula for the area of triangle T(t).
    """

    # Constants from the problem statement
    R = 10.0  # Radius of the circle and side length of the hexagon
    v = 1.0   # Speed of the triangle's vertices

    # --- Derivation of Coefficients ---
    # The area of the triangle T(t) is independent of the hexagon's rotation.
    # Due to symmetry, the triangle T(t) is always equilateral.
    # The area of an equilateral triangle with side length L is Area = (sqrt(3)/4) * L^2.
    # We need to find L^2 as a function of time t.

    # We set up a coordinate system with the hexagon's center at the origin.
    # The side length of the hexagon is s = R.
    # The displacement of each vertex of the triangle from the midpoint of its side
    # at time t is d = v * t = t.

    # The squared side length L^2 of the triangle can be derived from the geometry.
    # In a coordinate frame fixed to the hexagon, the vector connecting two
    # vertices of the triangle (e.g., one on side AB, one on side CD) has components:
    # Vx = -3 * R / 2
    # Vy = -sqrt(3) * R * t / s
    # Plugging in R=10 and s=10:
    # Vx = -15
    # Vy = -sqrt(3) * t

    # The squared side length L^2 is Vx^2 + Vy^2.
    # L_squared = (-15)^2 + (-sqrt(3) * t)^2
    # L_squared = 225 + 3 * t^2

    # From this derivation, we get the coefficients for the final area equation.
    # Area(t) = (sqrt(3)/4) * (225 + 3 * t^2)
    coeff_C1_num = "sqrt(3)"
    coeff_C1_den = 4
    coeff_C2 = 225.0
    coeff_C3 = 3.0

    # --- Final Result ---
    print("The area of the triangle T(t) is independent of the hexagon's rotation.")
    print("The triangle remains equilateral, and its area is a function of time t.")
    print("\nBased on the geometric derivation, the squared side length L(t)^2 is:")
    print(f"L(t)^2 = {int(coeff_C2)} + {int(coeff_C3)} * t^2")

    print("\nTherefore, the final equation for the area of the triangle T(t) is:")
    # The numbers in the final equation are C1_num, C1_den, C2, and C3.
    print(f"Area(t) = ({coeff_C1_num} / {coeff_C1_den}) * ({int(coeff_C2)} + {int(coeff_C3)} * t^2)")

solve_triangle_area()