import numpy as np

def solve():
    """
    Solves the problem by programmatically checking the Delaunay property for each triangulation.
    """
    # Estimated point coordinates based on the image, scaled for clarity.
    # Points P = {O1, O2, O3, O4, O5, I1, I2, I3}
    # Indices: O1=0, O2=1, O3=2, O4=3, O5=4, I1=5, I2=6, I3=7
    points = np.array([
        [0.0, 1.0],   # O1
        [0.9, 0.4],   # O2
        [0.6, -0.8],  # O3
        [-0.6, -0.8], # O4
        [-0.9, 0.4],  # O5
        [0.0, -0.2],  # I1
        [-0.5, 0.0],  # I2
        [0.5, 0.0],   # I3
    ])
    point_names = ["O1", "O2", "O3", "O4", "O5", "I1", "I2", "I3"]

    # Define the triangles for each figure.
    T_A = [
        (0, 1, 7), (0, 4, 6), (0, 6, 7), (5, 6, 7), (1, 2, 7),
        (4, 3, 6), (5, 2, 7), (5, 3, 6), (5, 2, 3)
    ]
    T_B = [
        (0, 1, 7), (0, 4, 6), (1, 2, 7), (2, 3, 5), (3, 4, 6),
        (0, 5, 6), (0, 5, 7), (5, 7, 2), (5, 6, 3)
    ]
    T_C = [
        (0, 1, 7), (0, 4, 6), (1, 2, 7), (3, 4, 6), (2, 3, 5),
        (0, 5, 7), (0, 5, 6), (5, 6, 3), (5, 7, 2)
    ]
    
    triangulations = {'A': T_A, 'B': T_B, 'C': T_C}
    violating_triangulations = []

    print("Analyzing triangulations for Delaunay property violations...")
    
    for name, tris in triangulations.items():
        is_delaunay = True
        for tri_indices in tris:
            p1_idx, p2_idx, p3_idx = tri_indices
            
            p1 = points[p1_idx]
            p2 = points[p2_idx]
            p3 = points[p3_idx]
            
            # Ensure CCW orientation for a positive determinant sign for "in"
            orient = (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[0] - p1[0])
            if orient < 0:
                p2_idx, p3_idx = p3_idx, p2_idx
                
            p1, p2, p3 = points[p1_idx], points[p2_idx], points[p3_idx]

            for p4_idx in range(len(points)):
                if p4_idx in [p1_idx, p2_idx, p3_idx]:
                    continue
                
                p4 = points[p4_idx]
                
                matrix = np.array([
                    [p1[0], p1[1], p1[0]**2 + p1[1]**2, 1],
                    [p2[0], p2[1], p2[0]**2 + p2[1]**2, 1],
                    [p3[0], p3[1], p3[0]**2 + p3[1]**2, 1],
                    [p4[0], p4[1], p4[0]**2 + p4[1]**2, 1]
                ])
                
                det = np.linalg.det(matrix)

                if det > 1e-9: # If determinant > 0, p4 is inside the circumcircle
                    print(f"\nViolation found in Triangulation {name}:")
                    tri_names = f"({point_names[p1_idx]}, {point_names[p2_idx]}, {point_names[p3_idx]})"
                    print(f"The circumcircle of triangle {tri_names} contains point {point_names[p4_idx]}.")
                    print("The in-circle test relies on the sign of a determinant.")
                    print("The equation is: det(M) > 0, where M is the 4x4 matrix of point coordinates.")
                    print("Matrix M:")
                    for row_idx, p_idx in enumerate([p1_idx, p2_idx, p3_idx, p4_idx]):
                        p = points[p_idx]
                        print(f"[{p[0]:6.2f} {p[1]:6.2f} {p[0]**2+p[1]**2:6.2f} {1:6.2f}]")
                    print(f"Calculated determinant: {det:.4f}")
                    print("Since the determinant is positive, the property is violated.")
                    is_delaunay = False
                    break
            if not is_delaunay:
                break
        
        if not is_delaunay:
            violating_triangulations.append(name)
            
    if not violating_triangulations:
        print("\nNo violations found in triangulations A, B, or C.")
    else:
        print(f"\nFinal Answer: The triangulation(s) violating the property are: {', '.join(violating_triangulations)}")

solve()