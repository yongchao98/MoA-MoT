import sympy

def solve_cone_distance():
    """
    Calculates the furthest distance one can travel from a point P on the base
    of a specific cone, traveling only on its surface.
    """

    # Define d as a positive symbolic variable, representing the diameter and slant height.
    d = sympy.Symbol('d', positive=True)

    # The problem can be solved by unrolling the cone into a semicircle of radius d.
    # We place the apex at the origin (0,0).
    # We choose the cut to be along the generator VP, so P is represented by two points
    # on the 2D plane: P1 at (d, 0) and P2 at (-d, 0).
    # The furthest point Q must be equidistant from P1 and P2, so it lies on the y-axis.
    # The point Q must also be within the semicircle (x^2 + y^2 <= d^2, y>=0).
    # This means Q has coordinates (0, y) where 0 <= y <= d.
    # We want to maximize the distance from P1=(d, 0) to Q=(0, y).
    # Distance^2 = (d - 0)^2 + (0 - y)^2 = d^2 + y^2.
    # This is maximized when y is maximized, so y = d.
    # The furthest point Q is at (0, d).

    # Now, we calculate this maximum distance.
    P1_x, P1_y = d, 0
    Q_x, Q_y = 0, d
    
    # Using the distance formula: sqrt((x2-x1)^2 + (y2-y1)^2)
    distance_squared = (P1_x - Q_x)**2 + (P1_y - Q_y)**2
    max_distance = sympy.sqrt(distance_squared)

    print("The furthest point Q from P is found by unrolling the cone into a semicircle.")
    print("The starting point P is at one end of the semicircle's diameter, say at (d, 0).")
    print("The furthest point Q is at the top of the semicircle's arc, at (0, d).")
    print("\nWe calculate the distance using the distance formula:")
    print(f"Distance = sqrt( (d - 0)^2 + (0 - d)^2 )")
    print(f"         = sqrt( {d**2} + {(-d)**2} )")
    print(f"         = sqrt( {d**2} + {d**2} )")
    print(f"         = sqrt( {2*d**2} )")
    print(f"\nThe final expression for the furthest distance is: {max_distance}")

solve_cone_distance()