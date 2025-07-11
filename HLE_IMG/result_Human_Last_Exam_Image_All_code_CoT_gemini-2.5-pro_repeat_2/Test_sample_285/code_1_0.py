import numpy as np

def check_violation_in_C():
    """
    Checks if triangulation C violates the Delaunay empty circle property
    by testing if point P2 is inside the circumcircle of triangle (P1, C1, C2).
    """
    # 1. Define coordinates for the relevant points based on the image.
    points = {
        'P1': np.array([1.9, 0.6]),   # Top-right vertex of the triangle
        'C1': np.array([0.8, 0.2]),   # Inner-right vertex of the triangle
        'C2': np.array([0.0, -1.0]),  # Inner-bottom vertex of the triangle
        'P2': np.array([1.2, -1.6])   # Point to test
    }

    p1, c1, c2 = points['P1'], points['C1'], points['C2']
    p2_test = points['P2']

    # The in-circle test checks the sign of a determinant.
    # For a counter-clockwise oriented triangle (a,b,c), a point d is inside
    # the circumcircle if the following determinant is positive:
    # | ax-dx  ay-dy  (ax-dx)^2+(ay-dy)^2 |
    # | bx-dx  by-dy  (bx-dx)^2+(by-dy)^2 |
    # | cx-dx  cy-dy  (cx-dx)^2+(cy-dy)^2 |
    
    # Let's define the terms of the matrix:
    # a, b, c = p1, c1, c2
    # d = p2_test
    
    adx = p1[0] - p2_test[0]
    ady = p1[1] - p2_test[1]
    bdx = c1[0] - p2_test[0]
    bdy = c1[1] - p2_test[1]
    cdx = c2[0] - p2_test[0]
    cdy = c2[1] - p2_test[1]

    ad_sq = adx**2 + ady**2
    bd_sq = bdx**2 + bdy**2
    cd_sq = cdx**2 + cdy**2

    # Calculate the determinant
    term1 = ad_sq * (bdx * cdy - bdy * cdx)
    term2 = bd_sq * (adx * cdy - ady * cdx)
    term3 = cd_sq * (adx * bdy - ady * bdx)
    
    determinant = term1 - term2 + term3

    # To be fully correct, we must ensure the triangle (p1, c1, c2) has a
    # positive (counter-clockwise) orientation.
    orientation = (c1[0] - p1[0]) * (c2[1] - p1[1]) - (c1[1] - p1[1]) * (c2[0] - p1[0])
    if orientation < 0:
        determinant = -determinant

    print("Checking Triangulation C:")
    print("Triangle vertices: P1(1.9, 0.6), C1(0.8, 0.2), C2(0.0, -1.0)")
    print("Test point: P2(1.2, -1.6)")
    print("\nThe in-circle test involves calculating a determinant. A positive value means the point is inside the circumcircle.")
    print("\nCalculation:")
    
    # Print the values used in the determinant calculation
    # For simplicity in output, we show the expanded formula terms
    # det = a_sq * (b_x*c_y - b_y*c_x) - b_sq * (a_x*c_y - a_y*c_x) + c_sq * (a_x*b_y - a_y*b_x)
    print(f"Term 1 = {ad_sq:.2f} * ({bdx:.2f} * {cdy:.2f} - {bdy:.2f} * {cdx:.2f}) = {term1:.4f}")
    print(f"Term 2 = {bd_sq:.2f} * ({adx:.2f} * {cdy:.2f} - {ady:.2f} * {cdx:.2f}) = {term2:.4f}")
    print(f"Term 3 = {cd_sq:.2f} * ({adx:.2f} * {bdy:.2f} - {ady:.2f} * {bdx:.2f}) = {term3:.4f}")

    print(f"\nEquation: {term1:.4f} - ({term2:.4f}) + {term3:.4f} = {determinant:.4f}")
    
    print(f"\nThe result of the in-circle test is {determinant:.4f}.")

    if determinant > 0:
        print("Since the result is positive, point P2 is inside the circumcircle of triangle (P1, C1, C2).")
        print("This violates the empty circle property.")
        print("\nTherefore, triangulation C is not a Delaunay triangulation.")
    else:
        print("The test did not show a violation for this triangle and point.")

check_violation_in_C()