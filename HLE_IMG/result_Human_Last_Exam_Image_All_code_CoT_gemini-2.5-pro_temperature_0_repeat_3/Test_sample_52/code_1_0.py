import math

def solve():
    """
    Calculates the period of tiles on a hexagonal grid for four different cases.
    """
    # 1. Establish a consistent coordinate system.
    # Side length s=2. C13 at origin.
    # This system is derived by enforcing that shared vertices are equidistant
    # from the centers of the hexagons they connect.
    s = 2.0
    sqrt3 = math.sqrt(3)

    # 2. Assign coordinates to key points.
    p = {}
    # Centers of the three main hexagons
    p[13] = (0.0, 0.0)
    p[31] = (3.0, sqrt3)
    p[23] = (3.0, -sqrt3)

    # Vertices and midpoints are derived from the geometry.
    # For example, vertex 8 is shared by H13, H31, and H_top.
    # Its coordinates can be found by solving for a point equidistant from the three centers.
    # This process yields the following coordinates for key vertices:
    p[8] = (1.0, sqrt3)
    p[9] = (-1.0, sqrt3)
    p[7] = (1.0, -sqrt3)
    p[12] = (-1.0, -sqrt3)

    # Midpoints lie on the center of edges connecting vertices.
    # Point 10 is the midpoint of the edge between vertex 8 and 9.
    p[10] = ((p[8][0] + p[9][0]) / 2, (p[8][1] + p[9][1]) / 2)
    # Point 4 is the midpoint of the edge between vertex 7 and 12.
    p[4] = ((p[7][0] + p[12][0]) / 2, (p[7][1] + p[12][1]) / 2)

    # For case 3, we need point 5. The diagram shows it's a vertex of H13
    # near points 4 and 6. The only vertex in this location is point 7.
    # So we assume the label '5' refers to the same geometric point as '7'.
    p[5] = p[7]

    # 3. & 4. Calculate the period for each case.

    print("Calculations:")

    # Case 1: 13, 31, 23
    p_start_1 = p[13]
    p_end_1 = p[23]
    v1 = (p_end_1[0] - p_start_1[0], p_end_1[1] - p_start_1[1])
    period1 = math.hypot(v1[0], v1[1])
    print("Case 1: Period = |p[23] - p[13]| = |({:.2f}, {:.2f}) - ({:.2f}, {:.2f})| = |({:.2f}, {:.2f})| = sqrt({:.2f}^2 + {:.2f}^2) = {:.4f}".format(
        p_end_1[0], p_end_1[1], p_start_1[0], p_start_1[1], v1[0], v1[1], v1[0], v1[1], period1))

    # Case 2: 10, 4, 23, 31
    p_start_2 = p[10]
    p_end_2 = p[31]
    v2 = (p_end_2[0] - p_start_2[0], p_end_2[1] - p_start_2[1])
    period2 = math.hypot(v2[0], v2[1])
    print("Case 2: Period = |p[31] - p[10]| = |({:.2f}, {:.2f}) - ({:.2f}, {:.2f})| = |({:.2f}, {:.2f})| = sqrt({:.2f}^2 + {:.2f}^2) = {:.4f}".format(
        p_end_2[0], p_end_2[1], p_start_2[0], p_start_2[1], v2[0], v2[1], v2[0], v2[1], period2))

    # Case 3: 5, 15, 17, 19, 21, 7
    p_start_3 = p[5]
    p_end_3 = p[7]
    # Since p[5] and p[7] are the same point, the path is a closed loop.
    # The period of a tiling made from a closed loop on a lattice is the
    # length of the shortest non-zero lattice vector.
    # The lattice vectors are T1=(3,sqrt3) and T2=(3,-sqrt3).
    # The shortest vector is T1-T2 = (0, 2*sqrt3).
    period3 = 2 * sqrt3
    print("Case 3: Path is a closed loop (p[5] == p[7]). Period is the shortest lattice vector = |(0, 2*sqrt(3))| = {:.4f}".format(period3))

    # Case 4: 4, 5, 14, 23, 18, 19, 21, 22, 31, 30, 9, 8, 13
    p_start_4 = p[4]
    p_end_4 = p[13]
    v4 = (p_end_4[0] - p_start_4[0], p_end_4[1] - p_start_4[1])
    period4 = math.hypot(v4[0], v4[1])
    print("Case 4: Period = |p[13] - p[4]| = |({:.2f}, {:.2f}) - ({:.2f}, {:.2f})| = |({:.2f}, {:.2f})| = sqrt({:.2f}^2 + {:.2f}^2) = {:.4f}".format(
        p_end_4[0], p_end_4[1], p_start_4[0], p_start_4[1], v4[0], v4[1], v4[0], v4[1], period4))

    # Final Answer
    final_answer = "{},{},{},{}".format(period1, period2, period3, period4)
    print("\nFinal Answer:")
    print(final_answer)
    return final_answer

final_answer_string = solve()
# The final answer needs to be wrapped in <<<>>>
# print(f"<<<{final_answer_string}>>>")