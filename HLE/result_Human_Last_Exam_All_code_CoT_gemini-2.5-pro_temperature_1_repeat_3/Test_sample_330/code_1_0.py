import numpy as np

def solve_square_reconstruction():
    """
    Finds the vertices of a square given one point on each side.
    """
    # Given points
    p1 = np.array([0.3511, 0.2027])
    p2 = np.array([0.6753, 0.8303])
    p3 = np.array([-0.2845, 0.9905])
    p4 = np.array([-0.128, 0.2218])

    # Step 1 & 2: Assume pairing (p1, p3) and (p2, p4) are on opposite sides.
    # Calculate vector v from p1 to p3
    v = p3 - p1
    
    # Rotate v by +90 degrees (counter-clockwise) to get vector p
    # Rotation matrix for +90 is [[0, -1], [1, 0]]
    # p = [-v_y, v_x]
    p = np.array([-v[1], v[0]])

    # Calculate auxiliary point X
    X = p2 + p

    # Step 3: Find the lines of the square
    # The first line (L_da) passes through p4 and X.
    # Its direction vector is d_vec.
    d_vec = X - p4
    d_unit = d_vec / np.linalg.norm(d_vec)
    
    # The perpendicular direction vector
    d_perp = np.array([-d_unit[1], d_unit[0]])

    # Now we find the vertices by intersecting the lines.
    # A vertex is the intersection of two perpendicular lines.
    # We can find the parameters for the parametric line equations.
    # Vertex A = intersection of line through p1 (dir d_perp) and line through p4 (dir d_unit)
    # p1 + u*d_perp = p4 + t*d_unit => u*d_perp - t*d_unit = p4 - p1
    # We can solve for u by taking dot product with d_perp
    u_A = np.dot(p4 - p1, d_perp)
    A = p1 + u_A * d_perp

    # Vertex B = intersection of line through p1 (dir d_perp) and line through p2 (dir d_unit)
    u_B = np.dot(p2 - p1, d_perp)
    B = p1 + u_B * d_perp
    
    # Vertex C = intersection of line through p3 (dir d_perp) and line through p2 (dir d_unit)
    u_C = np.dot(p2 - p3, d_perp)
    C = p3 + u_C * d_perp
    
    # Vertex D = intersection of line through p3 (dir d_perp) and line through p4 (dir d_unit)
    u_D = np.dot(p4 - p3, d_perp)
    D = p3 + u_D * d_perp

    vertices = [A, B, C, D]
    
    # Sort vertices by x-coordinate
    vertices.sort(key=lambda coord: coord[0])

    # Print the results
    print("The coordinates of the vertices of the square are:")
    for vertex in vertices:
        print(f"({vertex[0]:.2f}, {vertex[1]:.2f})")
    
    # Create the final answer string for extraction
    final_answer_str = ", ".join([f"({v[0]:.2f}, {v[1]:.2f})" for v in vertices])
    return final_answer_str

# Execute the function and capture the final answer string
final_answer = solve_square_reconstruction()

# The final answer in the required format for the system.
# The code block above already printed the human-readable output.
# This is for automated checking.
# print(f"\n<<<{final_answer}>>>")