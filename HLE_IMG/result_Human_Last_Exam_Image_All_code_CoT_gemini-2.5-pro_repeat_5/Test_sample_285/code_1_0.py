import math

def calculate_circumcircle_and_check_point(p1, p2, p3, p_test):
    """
    Calculates the circumcircle of a triangle (p1, p2, p3) and checks if p_test is inside.
    It prints the steps of the calculation.
    """
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = p3
    xp, yp = p_test

    print(f"Analyzing triangle with vertices P1{p1}, P2{p2}, P3{p3}.")
    print(f"Checking if test point P_test{p_test} is inside its circumcircle.\n")

    # The denominator of the circumcenter formulas.
    # If D is 0, the points are collinear.
    D = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
    if D == 0:
        print("The triangle vertices are collinear, no circumcircle exists.")
        return

    # Formulas for the circumcenter coordinates (Cx, Cy)
    Cx = (1/D) * ((x1**2 + y1**2) * (y2 - y3) + (x2**2 + y2**2) * (y3 - y1) + (x3**2 + y3**2) * (y1 - y2))
    Cy = (1/D) * ((x1**2 + y1**2) * (x3 - x2) + (x2**2 + y2**2) * (x1 - x3) + (x3**2 + y3**2) * (x2 - x1))

    # Calculate radius squared (R^2)
    R2 = (x1 - Cx)**2 + (y1 - Cy)**2

    # Calculate distance squared from the test point to the circumcenter (d^2)
    d2 = (xp - Cx)**2 + (yp - Cy)**2

    print(f"Equation for Circumcenter X-coordinate (Cx):")
    print(f"Cx = (1/{D:.2f}) * (({x1**2 + y1**2})*({y2 - y3}) + ({x2**2 + y2**2})*({y3 - y1}) + ({x3**2 + y3**2})*({y1 - y2}))")
    print(f"Cx = {Cx:.2f}\n")

    print(f"Equation for Circumcenter Y-coordinate (Cy):")
    print(f"Cy = (1/{D:.2f}) * (({x1**2 + y1**2})*({x3 - x2}) + ({x2**2 + y2**2})*({x1 - x3}) + ({x3**2 + y3**2})*({x2 - x1}))")
    print(f"Cy = {Cy:.2f}\n")

    print(f"Circumcenter is C = ({Cx:.2f}, {Cy:.2f})")

    print(f"Equation for Radius Squared (R^2):")
    print(f"R^2 = ({x1} - {Cx:.2f})^2 + ({y1} - {Cy:.2f})^2")
    print(f"R^2 = {R2:.2f}\n")

    print(f"Equation for Distance Squared from P_test to C (d^2):")
    print(f"d^2 = ({xp} - {Cx:.2f})^2 + ({yp} - {Cy:.2f})^2")
    print(f"d^2 = {d2:.2f}\n")

    print("Checking if the point is inside the circle (d^2 < R^2):")
    print(f"{d2:.2f} < {R2:.2f} ?")

    if d2 < R2:
        print("True. The point is inside the circumcircle.")
        print("\nConclusion: Triangulation B violates the Delaunay empty circle property.")
    else:
        print("False. The point is not inside the circumcircle.")

# Define coordinates for points from graph B based on visual estimation.
# P1: top-left, P2: top-right, P6: inner-right, P5: inner-left
P1 = (-10, 10)
P2 = (10, 10)
P6 = (4, 0)
P5 = (-3, 0) # The test point

# The triangle from graph B that we suspect violates the property
triangle_T = (P1, P2, P6)

# Run the check
calculate_circumcircle_and_check_point(triangle_T[0], triangle_T[1], triangle_T[2], P5)