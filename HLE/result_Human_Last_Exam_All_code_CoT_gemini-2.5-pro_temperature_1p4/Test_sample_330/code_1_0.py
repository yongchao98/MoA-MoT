import numpy as np

def solve_and_print_square_vertices():
    """
    Calculates and prints the vertices of a square given four points on its sides.
    """
    # Given points
    points = np.array([
        [0.3511, 0.2027],  # P1
        [0.6753, 0.8303],  # P2
        [-0.2845, 0.9905], # P3
        [-0.128, 0.2218]   # P4
    ])
    
    P1, P2, P3, P4 = points[0], points[1], points[2], points[3]

    # Assume P1/P3 are opposite, P2/P4 are opposite.
    # This leads to two possible orientations for the square. We choose the one
    # that results in a valid configuration (points on side segments).
    
    # Vectors between opposite points
    dx13 = P1[0] - P3[0]
    dy13 = P1[1] - P3[1]
    dx24 = P2[0] - P4[0]
    dy24 = P2[1] - P4[1]

    # Calculate the tangent of the angle 'theta' for the sides' orientation
    # m = tan(theta) = (dy13 + dx24) / (dx13 - dy24)
    m = (dy13 + dx24) / (dx13 - dy24)
    
    # Calculate sin(theta) and cos(theta)
    theta = np.arctan(m)
    c = np.cos(theta)  # cos(theta)
    s = np.sin(theta)  # sin(theta)

    # The four lines of the square are:
    # L_AB, passes through P1, direction (c, s) -> normal (-s, c)
    # L_BC, passes through P2, direction (-s, c) -> normal (c, s)
    # L_CD, passes through P3, direction (c, s) -> normal (-s, c)
    # L_DA, passes through P4, direction (-s, c) -> normal (c, s)
    # Line equation form: normal_x * x + normal_y * y = constant
    
    # Constants for the line equations
    # L_AB: -s*x + c*y = -s*P1x + c*P1y
    b_ab = -s * P1[0] + c * P1[1]
    # L_BC: c*x + s*y = c*P2x + s*P2y
    b_bc =  c * P2[0] + s * P2[1]
    # L_CD: -s*x + c*y = -s*P3x + c*P3y
    b_cd = -s * P3[0] + c * P3[1]
    # L_DA: c*x + s*y = c*P4x + s*P4y
    b_da =  c * P4[0] + s * P4[1]

    # Vertices are the intersections of these lines
    # Vertex A is the intersection of L_DA and L_AB
    # Vertex B is the intersection of L_AB and L_BC
    # Vertex C is the intersection of L_BC and L_CD
    # Vertex D is the intersection of L_CD and L_DA

    # We solve the 2x2 system for each vertex:
    # For A (L_DA and L_AB):
    #  c*x + s*y = b_da
    # -s*x + c*y = b_ab
    # Solution: x = c*b_da - s*b_ab, y = s*b_da + c*b_ab
    
    vA = np.array([c * b_da - s * b_ab, s * b_da + c * b_ab])
    vB = np.array([c * b_bc - s * b_ab, s * b_bc + c * b_ab])
    vC = np.array([c * b_bc - s * b_cd, s * b_bc + c * b_cd])
    vD = np.array([c * b_da - s * b_cd, s * b_da + c * b_cd])
    
    vertices = [vA, vB, vC, vD]
    
    print("The equations for the vertices are derived from the side lines of the square.")
    print("Below are the systems of linear equations solved for each vertex (A, B, C, D):\n")

    # Function to print the system of equations for a vertex
    def print_system(name, n1x, n1y, b1, n2x, n2y, b2):
        print(f"To find vertex {name}:")
        print(f"  {n1x:.4f}*x + {n1y:.4f}*y = {b1:.4f}")
        print(f"  {n2x:.4f}*x + {n2y:.4f}*y = {b2:.4f}\n")

    # Print equations for each vertex
    # A = L_DA intersection L_AB
    print_system('A', c, s, b_da, -s, c, b_ab)
    # B = L_AB intersection L_BC
    print_system('B', -s, c, b_ab, c, s, b_bc)
    # C = L_BC intersection L_CD
    print_system('C', c, s, b_bc, -s, c, b_cd)
    # D = L_CD intersection L_DA
    print_system('D', -s, c, b_cd, c, s, b_da)

    # Sort vertices by x-coordinate
    vertices.sort(key=lambda v: v[0])
    
    print("The coordinates of the vertices of the square, sorted by x-coordinate, are:")
    final_coords = []
    for v in vertices:
        final_coords.append(f"({v[0]:.2f},{v[1]:.2f})")
    print(", ".join(final_coords))

    return final_coords


# Execute the function to get the answer
final_coordinates = solve_and_print_square_vertices()
# Format for the final answer block
answer_string = ", ".join(final_coordinates)
print(f"\n<<<{answer_string}>>>")