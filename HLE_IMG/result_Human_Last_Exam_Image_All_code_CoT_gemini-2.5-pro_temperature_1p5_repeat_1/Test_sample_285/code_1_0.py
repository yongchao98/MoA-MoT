import numpy as np

def get_circumcircle_details(p1, p2, p3):
    """
    Calculates the circumcenter and circumradius squared of a triangle.
    Uses the formulas from Wikipedia: https://en.wikipedia.org/wiki/Circumscribed_circle
    """
    D = 2 * (p1[0] * (p2[1] - p3[1]) + p2[0] * (p3[1] - p1[1]) + p3[0] * (p1[1] - p2[1]))
    if D == 0:
        return None, None # Collinear points
        
    p1_sq = p1[0]**2 + p1[1]**2
    p2_sq = p2[0]**2 + p2[1]**2
    p3_sq = p3[0]**2 + p3[1]**2
    
    ux = (1 / D) * (p1_sq * (p2[1] - p3[1]) + p2_sq * (p3[1] - p1[1]) + p3_sq * (p1[1] - p2[1]))
    uy = (1 / D) * (p1_sq * (p3[0] - p2[0]) + p2_sq * (p1[0] - p3[0]) + p3_sq * (p2[0] - p1[0]))
    
    center = np.array([ux, uy])
    radius_sq = (p1[0] - ux)**2 + (p1[1] - uy)**2
    
    return center, radius_sq

def check_triangulation_b():
    """
    Checks if triangulation B violates the empty circle property.
    """
    # Estimated coordinates for the 7 points based on the image
    # P0=top, P1=top-right, P2=bottom-right, P3=bottom-left, P4=top-left, 
    # P5=inner-right, P6=inner-left
    points = {
        'P0': np.array([50, 95]),
        'P1': np.array([90, 65]),
        'P2': np.array([80, 10]),  # Point for triangle
        'P3': np.array([20, 10]),  # Point for triangle
        'P4': np.array([10, 65]),
        'P5': np.array([65, 45]),  # Point for triangle
        'P6': np.array([35, 45])   # Point to check
    }

    # In triangulation B, one of the triangles is (P3, P2, P5)
    p_triangle1 = points['P3']
    p_triangle2 = points['P2']
    p_triangle3 = points['P5']
    
    # We check if its circumcircle contains any other point.
    # The suspect is P6, as P2,P3,P5,P6 form an isosceles trapezoid.
    point_to_check = points['P6']

    # Calculate circumcircle for triangle (P3, P2, P5)
    center, radius_sq = get_circumcircle_details(p_triangle1, p_triangle2, p_triangle3)

    if center is None:
        print("The chosen triangle points are collinear.")
        return

    # Calculate distance squared from the point to check to the circumcenter
    dist_sq = np.sum((point_to_check - center)**2)

    print("Analysis for Triangulation B:")
    print("-" * 30)
    print(f"Triangle vertices: P3={p_triangle1}, P2={p_triangle2}, P5={p_triangle3}")
    print(f"Circumcenter: ({center[0]:.2f}, {center[1]:.2f})")
    print(f"Circumradius squared: {radius_sq:.2f}")
    print("-" * 30)
    print(f"Point to check: P6={point_to_check}")
    print(f"Distance squared from P6 to center: {dist_sq:.2f}")
    print("-" * 30)

    # A tiny tolerance is used to account for floating point inaccuracies.
    tolerance = 1e-9
    if dist_sq < radius_sq - tolerance:
        print("Result: P6 is INSIDE the circumcircle.")
        print("Conclusion: Triangulation B violates the empty circle property.")
    elif dist_sq <= radius_sq + tolerance:
        print("Result: P6 is ON the circumcircle (a degenerate case).")
        print("Conclusion: Triangulation B violates the empty circle property as the circle is not empty.")
    else:
        print("Result: P6 is OUTSIDE the circumcircle.")
        print("Conclusion: This specific check does not show a violation for B.")

check_triangulation_b()
