import numpy as np

def solve():
    """
    This script checks which triangulation violates the empty circle property of Delaunay triangulations.
    We will analyze triangulation B.
    """
    # Step 1: Define the points based on visual estimation from the image.
    # The exact values are not critical, their relative positions are.
    points = {
        'P0': np.array([0, 0]),
        'P1': np.array([20, 0]),
        'P2': np.array([25, 20]),
        'P3': np.array([10, 25]),
        'P4': np.array([-5, 20]),
        'P5': np.array([2, 10]),
        'P6': np.array([10, 5]),
        'P7': np.array([18, 10]),
    }

    # Step 2: In Graph B, identify a triangle to test. Let's choose T(P3, P5, P6).
    # We will check if its circumcircle contains point P7.
    triangle_vertices_names = ['P3', 'P5', 'P6']
    triangle_vertices_coords = [points[name] for name in triangle_vertices_names]
    point_to_check_name = 'P7'
    point_to_check_coord = points[point_to_check_name]
    
    p1, p2, p3 = triangle_vertices_coords

    # Step 3: Calculate the circumcenter and squared radius of the triangle.
    # Using the formula based on the intersection of perpendicular bisectors.
    D = 2 * (p1[0] * (p2[1] - p3[1]) + p2[0] * (p3[1] - p1[1]) + p3[0] * (p1[1] - p2[1]))
    
    if D == 0:
        print("The chosen vertices are collinear.")
        return

    p1_sq = np.sum(p1**2)
    p2_sq = np.sum(p2**2)
    p3_sq = np.sum(p3**2)
    
    ux = (p1_sq * (p2[1] - p3[1]) + p2_sq * (p3[1] - p1[1]) + p3_sq * (p1[1] - p2[1])) / D
    uy = (p1_sq * (p3[0] - p2[0]) + p2_sq * (p1[0] - p3[0]) + p3_sq * (p2[0] - p1[0])) / D
    
    center = np.array([ux, uy])
    radius_sq = np.sum((p1 - center)**2)

    # Step 4: Check if the point is inside the circle.
    dist_sq = np.sum((point_to_check_coord - center)**2)

    # Step 5: Print the detailed results of the check.
    print("Analyzing triangulation B:")
    print(f"Chosen triangle vertices: {triangle_vertices_names[0]}{p1}, {triangle_vertices_names[1]}{p2}, {triangle_vertices_names[2]}{p3}")
    print(f"Point to check: {point_to_check_name}{point_to_check_coord}")
    print("-" * 30)
    print(f"Circumcenter of the triangle: ({center[0]:.2f}, {center[1]:.2f})")
    print(f"Squared radius of the circumcircle: {radius_sq:.2f}")
    print(f"Squared distance from {point_to_check_name} to the circumcenter: {dist_sq:.2f}")
    print("-" * 30)
    
    # Final conclusion based on the "empty circle" property test
    print("The empty circle property states that no point should be inside a triangle's circumcircle.")
    print(f"Test: Is squared_distance < squared_radius?")
    print(f"Equation: Is {dist_sq:.2f} < {radius_sq:.2f}?")
    
    if dist_sq < radius_sq:
        print("Result: Yes, the condition is true.")
        print(f"Conclusion: Point {point_to_check_name} lies inside the circumcircle of triangle ({', '.join(triangle_vertices_names)}).")
        print("Therefore, triangulation B violates the empty circle property.")
    else:
        print("Result: No, the condition is false.")
        print("This specific check does not show a violation.")

solve()