import numpy as np

def get_circumcircle_and_test_point():
    """
    Calculates the circumcircle for a specific triangle from Triangulation C and 
    checks if it contains another specific point, demonstrating a violation.
    """
    
    # Define points that capture the visual properties of triangulation C.
    # The triangle (P_top, P_tr, P_ir) is made to be very "skinny".
    P_top = (0, 100) 
    P_tr  = (9, 11)   # Placed very close to the line segment between P_top and P_ir
    P_ir  = (10, 0)
    P_il  = (-15, 2)  # A point to test for inclusion
    
    # The triangle in question from triangulation C
    triangle_vertices = [P_top, P_tr, P_ir]
    p1, p2, p3 = triangle_vertices

    print(f"Checking the circumcircle of triangle with vertices P_top={p1}, P_tr={p2}, P_ir={p3}")

    # Calculate circumcenter (cx, cy) and radius squared (r_sq)
    # Using the determinant formula for circumcircle calculation.
    D = 2 * (p1[0] * (p2[1] - p3[1]) + p2[0] * (p3[1] - p1[1]) + p3[0] * (p1[1] - p2[1]))

    # Avoid division by zero for collinear points
    if abs(D) < 1e-9:
        print("The points of the triangle are collinear.")
        return

    p1_sq = p1[0]**2 + p1[1]**2
    p2_sq = p2[0]**2 + p2[1]**2
    p3_sq = p3[0]**2 + p3[1]**2
    
    cx = (p1_sq * (p2[1] - p3[1]) + p2_sq * (p3[1] - p1[1]) + p3_sq * (p1[1] - p2[1])) / D
    cy = (p1_sq * (p3[0] - p2[0]) + p2_sq * (p1[0] - p3[0]) + p3_sq * (p2[0] - p1[0])) / D
    
    r_sq = (p1[0] - cx)**2 + (p1[1] - cy)**2

    print(f"Calculated Circumcenter: ({cx:.2f}, {cy:.2f})")
    print(f"Calculated Radius Squared: {r_sq:.2f}")

    # Check if the point P_il is inside this circumcircle.
    point_to_test = P_il
    print(f"\nTesting if point P_il={point_to_test} is inside the circumcircle...")

    dist_sq = (point_to_test[0] - cx)**2 + (point_to_test[1] - cy)**2
    
    print(f"Squared distance from center to P_il: {dist_sq:.2f}")
    
    if dist_sq < r_sq - 1e-9: # Using a small epsilon for float comparison
        print("\nResult: The point P_il is INSIDE the circumcircle.")
        print("Therefore, triangulation C violates the empty circle property of Delaunay triangulations.")
    else:
        print("\nResult: The point P_il is NOT inside the circumcircle.")

get_circumcircle_and_test_point()