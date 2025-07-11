import numpy as np

def solve_delaunay_violation():
    """
    Analyzes triangulation C to check for violations of the empty circle property.
    """
    # Step 1: Define the coordinates of the points based on the layout in Diagram C.
    points = {
        'Pt': np.array([0, 10]),      # Top point
        'Ptl': np.array([-2, 8]),     # Top-left point
        'Ptr': np.array([2, 8]),      # Top-right point
        'Pbl': np.array([-8, -8]),    # Bottom-left point
        'Pbr': np.array([8, -8]),     # Bottom-right point
        'I_high': np.array([0, 2]),   # Upper interior point
        'I_low': np.array([0, -2]),   # Lower interior point
    }

    # Step 2: Select a suspicious triangle from triangulation C and a point to check.
    # The triangle formed by (Pt, Ptl, I_high) is very skinny and a good candidate.
    # We will check if its circumcircle contains the point Ptr.
    p1, p2, p3 = points['Pt'], points['Ptl'], points['I_high']
    point_to_check = points['Ptr']
    
    print("Analyzing Triangulation C for Delaunay Property Violation")
    print("="*60)
    print("The empty circle property states that the circumcircle of any triangle in a Delaunay triangulation must contain no other points from the set.")
    print("\nWe will test the triangle formed by points Pt, Ptl, and I_high.")
    print(f"Vertices of the triangle: Pt={p1}, Ptl={p2}, I_high={p3}")

    # Step 3: Calculate the circumcenter and radius squared of the triangle.
    # We solve the system of linear equations for the perpendicular bisectors.
    # 2(x2-x1)x + 2(y2-y1)y = x2^2+y2^2 - x1^2-y1^2
    # 2(x3-x2)x + 2(y3-y2)y = x3^2+y3^2 - x2^2-y2^2
    
    A = np.array([
        [2 * (p2[0] - p1[0]), 2 * (p2[1] - p1[1])],
        [2 * (p3[0] - p2[0]), 2 * (p3[1] - p2[1])]
    ])
    
    b = np.array([
        p2[0]**2 + p2[1]**2 - p1[0]**2 - p1[1]**2,
        p3[0]**2 + p3[1]**2 - p2[0]**2 - p2[1]**2
    ])
    
    try:
        center = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        print("\nThe points are collinear and do not form a triangle.")
        return

    radius_sq = (p1[0] - center[0])**2 + (p1[1] - center[1])**2

    print(f"\nThe calculated circumcenter of the triangle is ({center[0]:.1f}, {center[1]:.1f}).")
    print(f"The radius squared of the circumcircle is {radius_sq:.1f}.")

    # Step 4: Check if the point_to_check lies inside the circle.
    dist_sq = (point_to_check[0] - center[0])**2 + (point_to_check[1] - center[1])**2
    
    print(f"\nNow, let's check if the point Ptr={point_to_check} is inside this circle.")
    print(f"The equation of the circumcircle is (x - {center[0]:.1f})^2 + (y - {center[1]:.1f})^2 = {radius_sq:.1f}")
    
    print(f"To check Ptr, we calculate its distance squared from the center:")
    print(f"  d^2 = ({point_to_check[0]} - {center[0]:.1f})^2 + ({point_to_check[1]} - {center[1]:.1f})^2 = {dist_sq:.1f}")

    print(f"\nComparing the distance squared to the radius squared: {dist_sq:.1f} < {radius_sq:.1f}")

    if dist_sq < radius_sq:
        print("\nConclusion: The point Ptr lies INSIDE the circumcircle.")
        print("Therefore, triangulation C violates the empty circle property of Delaunay triangulations.")
    else:
        print("\nConclusion: The point Ptr lies outside or on the circumcircle.")
        print("This specific triangle does not violate the property with this point.")

solve_delaunay_violation()