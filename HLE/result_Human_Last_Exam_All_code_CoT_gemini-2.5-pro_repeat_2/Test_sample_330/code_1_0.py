import numpy as np

def solve_square_from_points():
    """
    Finds the vertices of a square given one point on each of its four sides.
    """
    # 1. Define the four points in cyclic order
    P1 = np.array([0.3511, 0.2027])
    P2 = np.array([0.6753, 0.8303])
    P3 = np.array([-0.2845, 0.9905])
    P4 = np.array([-0.128, 0.2218])

    # 2. Calculate vectors between opposite points
    v1 = P3 - P1
    v2 = P4 - P2

    # Define a 90-degree counter-clockwise rotation matrix function
    def R90(v):
        return np.array([-v[1], v[0]])

    # 3. Determine the two possible direction vectors for the sides
    # There are two solutions, one is typically the correct one and the other is degenerate.
    # The side direction vectors are parallel to v1 ± R90(v2).
    dir_a = v1 + R90(v2)
    dir_b = v1 - R90(v2)

    # To choose the correct one, we can compare the resulting square side lengths.
    # The side length is proportional to |v1 . v2| / |dir|.
    # We choose the direction with the smaller magnitude, which gives a larger side length.
    if np.linalg.norm(dir_a) < np.linalg.norm(dir_b):
        side_dir = dir_a
    else:
        side_dir = dir_b

    # 4. Define the lines for the sides of the square
    # Normal vector to side AB (and CD)
    n_AC = R90(side_dir)
    # Normal vector to side BC (and DA)
    n_BD = -side_dir

    # Calculate the constants c for the line equations n.x = c
    c_A = np.dot(n_AC, P1) # Line AB
    c_C = np.dot(n_AC, P3) # Line CD
    c_B = np.dot(n_BD, P2) # Line BC
    c_D = np.dot(n_BD, P4) # Line DA

    # 5. Calculate the vertices by solving systems of linear equations
    # The system for intersection of ax+by=c1 and -bx+ay=c2 is M*v=r
    # where M = [[a,b],[-b,a]]
    M = np.array([n_AC, n_BD])

    # Vertex A is the intersection of Line DA and Line AB
    # n_BD.x = c_D and n_AC.x = c_A
    rhs_A = np.array([c_A, c_D])
    vertex_A = np.linalg.solve(M, rhs_A)

    # Vertex B is the intersection of Line AB and Line BC
    # n_AC.x = c_A and n_BD.x = c_B
    rhs_B = np.array([c_A, c_B])
    vertex_B = np.linalg.solve(M, rhs_B)

    # Vertex C is the intersection of Line BC and Line CD
    # n_BD.x = c_B and n_AC.x = c_C
    rhs_C = np.array([c_C, c_B])
    vertex_C = np.linalg.solve(M, rhs_C)

    # Vertex D is the intersection of Line CD and Line DA
    # n_AC.x = c_C and n_BD.x = c_D
    rhs_D = np.array([c_C, c_D])
    vertex_D = np.linalg.solve(M, rhs_D)

    vertices = [vertex_A, vertex_B, vertex_C, vertex_D]

    # 6. Sort vertices by x-coordinate
    vertices.sort(key=lambda v: v[0])

    # Print the results
    print("The coordinates of the vertices of the square are:")
    for v in vertices:
        print(f"({v[0]:.2f},{v[1]:.2f})")

solve_square_from_points()