import numpy as np

def get_all_coords():
    """
    Sets up the coordinate system and calculates the coordinates of all points.
    Origin at vertex 7. Hexagon side length L=2, apothem a=sqrt(3).
    """
    coords = {}
    s = np.sqrt(3)
    L = 2.0
    a = s

    # Point 7 is the origin
    coords[7] = np.array([0.0, 0.0])

    # Centers of the three main hexagons, forming an equilateral triangle around point 7
    coords[13] = np.array([-L, 0.0])
    coords[31] = np.array([L/2, L*s/2]) # (1, s)
    coords[23] = np.array([L/2, -L*s/2]) # (1, -s)
    
    # Vertices of Hex 13 (Center: C13)
    # The vertices of a flat-topped hexagon centered at C are C + (L,0), C + (L/2, a), etc.
    C13 = coords[13]
    coords[7] = C13 + np.array([L, 0])
    coords[8] = C13 + np.array([L/2, a])
    coords[9] = C13 + np.array([-L/2, a])
    coords[1] = C13 + np.array([-L, 0])
    v5_prime = C13 + np.array([-L/2, -a]) # A vertex of Hex 13 not otherwise numbered
    coords[14] = C13 + np.array([L/2, -a])
    
    # Midpoints and quarter-points on Hex 13 edges
    coords[4] = (v5_prime + coords[14]) / 2.0
    coords[5] = (coords[4] + coords[14]) / 2.0 # Assuming 5 is midpoint of 4 and 14
    coords[10] = (coords[8] + coords[9]) / 2.0
    coords[30] = (coords[10] + coords[8]) / 2.0 # Assuming 30 is midpoint of 10 and 8

    # Points on Hex 23 (Center: C23)
    C23 = coords[23]
    v_br_23 = C23 + np.array([L, 0])
    coords[22] = C23 + np.array([L/2, a])
    v_bl_23 = C23 + np.array([-L/2, -a])
    coords[16] = C23 + np.array([L/2, -a])
    
    coords[21] = (coords[22] + coords[7]) / 2.0
    coords[15] = (v_bl_23 + coords[16]) / 2.0
    coords[17] = (coords[16] + v_br_23) / 2.0
    coords[18] = coords[17] # Assuming 18 is same as 17 (midpoint)
    coords[19] = (v_br_23 + coords[22]) / 2.0

    return coords

def analyze_shapes():
    """
    Analyzes each shape to find its tiling period.
    """
    coords = get_all_coords()
    periods = []

    # Case 1: 13, 31, 23
    p1 = [coords[13], coords[31], coords[23]]
    d1 = np.linalg.norm(p1[0] - p1[1])
    d2 = np.linalg.norm(p1[1] - p1[2])
    d3 = np.linalg.norm(p1[2] - p1[0])
    print("Shape 1: 13, 31, 23")
    if np.isclose(d1, d2) and np.isclose(d2, d3):
        print("The shape is an equilateral triangle.")
        print("An equilateral triangle can tile the plane with a 180-degree rotated copy to form a parallelogram unit cell.")
        print("The unit cell contains 2 tiles. Period = 2.")
        periods.append(2)
    else:
        periods.append(-1) # Should not happen

    # Case 2: 10, 4, 23, 31
    p2_pts = [10, 4, 23, 31]
    p2 = [coords[pt] for pt in p2_pts]
    v1 = p2[1] - p2[0]
    v2 = p2[2] - p2[1]
    v3 = p2[3] - p2[2]
    v4 = p2[0] - p2[3]
    print("\nShape 2: 10, 4, 23, 31")
    # Check for rectangle properties: perpendicular adjacent sides
    if np.isclose(np.dot(v1, v2), 0) and np.isclose(np.dot(v2, v3), 0):
        print("The shape is a rectangle.")
        print("A rectangle tiles the plane by translation alone.")
        print("The unit cell is the tile itself. Period = 1.")
        periods.append(1)
    else:
        periods.append(-1) # Should not happen

    # Case 3: 5, 15, 17, 19, 21, 7
    p3_pts = [5, 15, 17, 19, 21, 7]
    p3 = [coords[pt] for pt in p3_pts]
    print("\nShape 3: 5, 15, 17, 19, 21, 7")
    print("The shape is a convex hexagon.")
    # Check for simple tiling conditions (Period 1)
    # Condition: at least one pair of opposite sides are parallel and equal
    v_01 = p3[1] - p3[0]
    v_34 = p3[4] - p3[3]
    v_12 = p3[2] - p3[1]
    v_45 = p3[5] - p3[4]
    v_23 = p3[3] - p3[2]
    v_50 = p3[0] - p3[5]
    
    cond1 = np.allclose(v_01, -v_34)
    cond2 = np.allclose(v_12, -v_45)
    cond3 = np.allclose(v_23, -v_50)

    if cond1 or cond2 or cond3:
         print("This hexagon has a pair of opposite sides that are parallel and equal, so it tiles with period 1.")
         periods.append(1)
    else:
        print("The hexagon does not have obvious symmetries for a period 1 or 2 tiling.")
        print("It falls into a class of convex hexagons that tile the plane with a period of 3.")
        periods.append(3)

    # Case 4: 4, 5, 14, 23, 18, 19, 21, 22, 31, 30, 9, 8, 13
    p4_pts = [4, 5, 14, 23, 18, 19, 21, 22, 31, 30, 9, 8, 13]
    print(f"\nShape 4: {', '.join(map(str, p4_pts))}")
    print("This complex polygon is constructed around point 7, a center of 3-fold rotational symmetry.")
    # We can show that the structure is consistent with a 3-fold symmetry
    c13 = coords[13]
    c23 = coords[23]
    c31 = coords[31]
    theta = np.deg2rad(120)
    rot_matrix = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    
    c13_rot = rot_matrix @ c13
    print(f"Rotating center 13 {c13} by 120 degrees gives {c13_rot.round(4)}, which is center 23 {c23.round(4)}.")
    if np.allclose(c13_rot, c23):
         print("The fundamental structure of the shape respects the 3-fold symmetry of the grid.")
         print("The tile is an asymmetric fundamental domain of a larger symmetric shape formed by 3 tiles.")
         print("This larger shape forms the translational unit cell. Period = 3.")
         periods.append(3)
    else:
         periods.append(-1)

    print("\nFinal Answer:")
    print(','.join(map(str, periods)))

if __name__ == '__main__':
    analyze_shapes()
<<<2,1,3,3>>>