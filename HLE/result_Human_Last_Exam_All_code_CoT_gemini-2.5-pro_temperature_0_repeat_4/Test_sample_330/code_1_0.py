import numpy as np

def solve_square_vertices():
    """
    Calculates the vertices of a square given four points on its sides.
    """
    # Let the four points be A, B, C, D, assuming they are on consecutive sides.
    A = np.array([0.3511, 0.2027])
    B = np.array([0.6753, 0.8303])
    C = np.array([-0.2845, 0.9905])
    D = np.array([-0.128, 0.2218])

    # There are two solutions. We'll find one by rotating B around the midpoint of AC by -90 degrees.
    # A +90 degree rotation would give the other solution.
    
    # 1. Find the midpoint of the segment AC
    M_ac = (A + C) / 2

    # 2. Find the vector from M to B
    vec_MB = B - M_ac

    # 3. Rotate the vector MB by -90 degrees.
    # A -90 degree rotation matrix is [[0, 1], [-1, 0]]
    vec_MB_rot = np.array([vec_MB[1], -vec_MB[0]])

    # 4. Find the rotated point B_star
    B_star = M_ac + vec_MB_rot

    # 5. The line passing through D and B_star is one side of the square.
    # Its direction vector determines the orientation of the square.
    d1 = B_star - D
    d1 = d1 / np.linalg.norm(d1) # Normalize the direction vector

    # 6. The other direction is perpendicular to d1.
    d2 = np.array([-d1[1], d1[0]])

    # 7. Define the four lines that form the square's sides.
    # A point and a direction vector define each line.
    # Line L_A passes through A, parallel to d2
    # Line L_B passes through B, parallel to d1
    # Line L_C passes through C, parallel to d2
    # Line L_D passes through D, parallel to d1

    def intersect(p1, v1, p2, v2):
        """
        Finds the intersection of two lines defined by a point and a direction vector.
        Line 1: p1 + t*v1
        Line 2: p2 + s*v2
        Solves the system: t*v1 - s*v2 = p2 - p1
        """
        A = np.array([[v1[0], -v2[0]], [v1[1], -v2[1]]])
        b = p2 - p1
        try:
            # Solve for [t, s]
            params = np.linalg.solve(A, b)
            # Return the intersection point using parameter t
            return p1 + params[0] * v1
        except np.linalg.LinAlgError:
            return None

    # 8. Calculate the vertices by finding the intersections of these lines.
    V4 = intersect(D, d1, C, d2)
    V3 = intersect(C, d2, B, d1)
    V2 = intersect(B, d1, A, d2)
    V1 = intersect(A, d2, D, d1)

    vertices = [V1, V2, V3, V4]

    # 9. Sort vertices by x-coordinate
    vertices.sort(key=lambda v: v[0])

    # 10. Print the results
    print("The coordinates of the vertices of the square are:")
    for v in vertices:
        print(f"({v[0]:.2f}, {v[1]:.2f})")

solve_square_vertices()
<<<
(-0.59, 0.54)
(-0.06, 1.33)
(0.21, -0.01)
(0.74, 0.78)
>>>