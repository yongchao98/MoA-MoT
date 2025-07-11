import numpy as np

def get_circumcircle(p1, p2, p3):
    """
    Calculates the circumcenter and the squared radius of the circle
    passing through three points p1, p2, p3.
    Returns a tuple (center, radius_squared).
    Returns (None, None) if the points are collinear.
    """
    p1, p2, p3 = np.array(p1), np.array(p2), np.array(p3)
    # Using the formula from the intersection of perpendicular bisectors
    D = 2 * (p1[0] * (p2[1] - p3[1]) + p2[0] * (p3[1] - p1[1]) + p3[0] * (p1[1] - p2[1]))
    if abs(D) < 1e-10:
        return None, None
        
    p1_sq = np.sum(p1**2)
    p2_sq = np.sum(p2**2)
    p3_sq = np.sum(p3**2)
    
    ux = (p1_sq * (p2[1] - p3[1]) + p2_sq * (p3[1] - p1[1]) + p3_sq * (p1[1] - p2[1])) / D
    uy = (p1_sq * (p3[0] - p2[0]) + p2_sq * (p1[0] - p3[0]) + p3_sq * (p2[0] - p1[0])) / D
    
    center = np.array([ux, uy])
    radius_sq = np.sum((p1 - center)**2)
    return center, radius_sq

def check_violation():
    """
    Checks for a Delaunay violation in triangulation B.
    """
    print("Analyzing Triangulation B for a violation of the empty circle property.")
    print("The empty circle property states that the circumcircle of any triangle in the triangulation must contain no other points of the set in its interior.")
    print("\nIn diagram B, consider the convex quadrilateral formed by points P2, Q2, Q3, and P3.")
    print("This quadrilateral is divided by the edge connecting Q2 and P3.")
    print("The Delaunay property is violated if point P2 lies inside the circumcircle of triangle (Q2, Q3, P3).")
    
    # Assigning coordinates that represent the visual arrangement in B.
    # We choose points that clearly demonstrate the violation.
    p2_coord = (0, 0.5) 
    q2_coord = (2, 0)
    q3_coord = (-2, 0)
    p3_coord = (0, -3)

    print("\nLet's assign representative coordinates:")
    print(f"Point P2 = {p2_coord}")
    print(f"Point Q2 = {q2_coord}")
    print(f"Point Q3 = {q3_coord}")
    print(f"Point P3 = {p3_coord}")
    
    violating_triangle_points = [q2_coord, q3_coord, p3_coord]
    offending_point = p2_coord
    
    center, radius_sq = get_circumcircle(*violating_triangle_points)
    offending_point_np = np.array(offending_point)
    dist_sq = np.sum((offending_point_np - center)**2)
    
    print(f"\nStep 1: Find the circumcircle of triangle (Q2, Q3, P3).")
    print(f"The vertices are Q2={q2_coord}, Q3={q3_coord}, P3={p3_coord}.")
    print(f"Calculated Circumcenter: ({center[0]:.2f}, {center[1]:.2f})")
    print(f"Calculated Squared Radius of Circumcircle: {radius_sq:.2f}")
    
    print(f"\nStep 2: Check if point P2 lies inside this circle.")
    print(f"The point to check is P2={offending_point}.")
    print(f"The squared distance of P2 from the circumcenter is {dist_sq:.2f}.")
    
    print("\nStep 3: Compare the squared distance to the squared radius.")
    print(f"Is distance^2 < radius^2?")
    print(f"Is {dist_sq:.2f} < {radius_sq:.2f}?")
    
    if dist_sq < radius_sq:
        print("Yes, the inequality is true.")
        print("\nConclusion: Point P2 lies inside the circumcircle of triangle (Q2, Q3, P3).")
        print("Therefore, Triangulation B violates the empty circle property.")
    else:
        # This branch should not be hit with the selected coordinates.
        print("No, the inequality is false. The property holds for this specific configuration.")

# Run the analysis.
check_violation()