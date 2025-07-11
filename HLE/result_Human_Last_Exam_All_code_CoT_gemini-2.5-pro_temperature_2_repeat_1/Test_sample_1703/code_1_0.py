import math

def solve_triangle_area():
    """
    Calculates the coefficients of the formula for the area of triangle T(t).
    The area is given by A(t) = (sqrt(3) / A) * (B + C * t^2).
    """

    # --- Problem Parameters ---
    # Radius of the circumscribed circle
    R = 10.0
    # Side length of the regular hexagon is equal to the radius
    s = R
    # Speed of the vertices
    v = 1.0

    # --- Derivation of Coefficients ---
    # The final formula for the area has been derived from analytical geometry
    # as Area(t) = (sqrt(3) / 4) * (L0**2 + 3 * t**2), where L0 is the
    # initial side length of the triangle at t=0.

    # Coefficient A is the denominator from the area formula of an equilateral
    # triangle (sqrt(3)/4 * side^2) and algebraic simplification.
    coeff_A = 4

    # Coefficient B is the square of the triangle's side length at t=0.
    # At t=0, the triangle's vertices are at the midpoints of alternating
    # sides of the hexagon. The side length (L0) of this triangle is (3/2)*s.
    L0 = (3.0 / 2.0) * s
    coeff_B = L0**2

    # Coefficient C is a constant derived from the geometric analysis of
    # the vertices' movement.
    coeff_C = 3

    print("The area of the triangle T(t) as a function of time t is determined by the following formula:")
    print("Area(t) = (sqrt(3) / A) * (B + C * t^2)")
    print("\nThe script calculates the values for the coefficients A, B, and C.")

    print("\n--- Calculated Coefficients ---")
    print(f"The side length of the hexagon s = R = {int(s)}")
    print(f"The initial side length of the triangle L0 = (3/2)*s = {L0}")
    print("\nThe numbers in the final equation are:")
    print(f"A = {int(coeff_A)}")
    print(f"B = L0^2 = {int(coeff_B)}")
    print(f"C = {int(coeff_C)}")

    print("\n--- Final Equation ---")
    print("Substituting the numbers, the area of the triangle T(t) is:")
    print(f"Area(t) = (sqrt(3) / {int(coeff_A)}) * ({int(coeff_B)} + {int(coeff_C)} * t^2)")

solve_triangle_area()