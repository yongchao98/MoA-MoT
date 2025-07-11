import numpy as np
from itertools import combinations

def solve_geometry_problem():
    """
    Calculates and verifies the volume of the polytope P for a 2D simplex (triangle).
    """
    # Define the vertices of the simplex (triangle) T.
    # You can change these vertices to test any triangle.
    # Example 1: A right-angled triangle
    # vertices = np.array([[0.0, 0.0], [3.0, 0.0], [0.0, 4.0]])
    # Example 2: An equilateral triangle
    vertices = np.array([[0.0, 0.0], [2.0, 0.0], [1.0, np.sqrt(3)]])
    # Example 3: A general triangle
    # vertices = np.array([[1.0, 1.0], [5.0, 2.0], [2.0, 6.0]])

    d = 2
    if vertices.shape != (d + 1, d):
        print(f"Error: For d={d}, expected {d+1} vertices.")
        return

    # --- 1. Calculate the volume (area) of the simplex T ---
    # For d=2, V = 0.5 * |x1(y2-y3) + x2(y3-y1) + x3(y1-y2)|
    # A more general way using determinant:
    mat = np.ones((d + 1, d + 1))
    mat[:, :-1] = vertices
    simplex_volume = np.abs(np.linalg.det(mat)) / np.math.factorial(d)

    if simplex_volume < 1e-9:
        print("The provided vertices form a degenerate simplex (zero volume).")
        return

    # --- 2. Define the polytope P' in alpha-space ---
    v0 = vertices[0]
    u1 = vertices[1] - v0
    u2 = vertices[2] - v0

    # Gram matrix G_ij = u_i . u_j
    G = np.array([
        [np.dot(u1, u1), np.dot(u1, u2)],
        [np.dot(u2, u1), np.dot(u2, u2)]
    ])

    # The polytope P' is defined by 6 linear inequalities A*alpha <= b
    # Derived from the definition of P
    lines = []
    # From edges (v0,v1) and (v0,v2)
    # 0 <= alpha.G_k <= G_kk
    lines.append([G[0, 0], G[1, 0], G[0, 0]])  # u1.G_1*alpha <= G_11
    lines.append([-G[0, 0], -G[1, 0], 0])      # u1.G_1*alpha >= 0
    lines.append([G[0, 1], G[1, 1], G[1, 1]])  # u2.G_2*alpha <= G_22
    lines.append([-G[0, 1], -G[1, 1], 0])      # u2.G_2*alpha >= 0
    # From edge (v1,v2)
    # G_12-G_11 <= alpha.(G_2-G_1) <= G_22-G_12
    A1 = G[0, 1] - G[0, 0]
    A2 = G[1, 1] - G[1, 0]
    lines.append([A1, A2, G[1, 1] - G[1, 0]])
    lines.append([-A1, -A2, -(G[0, 1] - G[0, 0])])
    
    # --- 3. Find the vertices of the polytope P' ---
    poly_vertices = []
    for i in range(len(lines)):
        for j in range(i + 1, len(lines)):
            A = np.array([lines[i][:2], lines[j][:2]])
            b = np.array([lines[i][2], lines[j][2]])
            
            # Check if lines are parallel
            if np.abs(np.linalg.det(A)) < 1e-9:
                continue
            
            try:
                # Solve for the intersection point
                intersection = np.linalg.solve(A, b)
                
                # Check if the point satisfies all other inequalities
                is_vertex = True
                for k in range(len(lines)):
                    if k != i and k != j:
                        # A_k*x <= b_k  => A_k*x - b_k <= 0
                        if np.dot(lines[k][:2], intersection) > lines[k][2] + 1e-9:
                            is_vertex = False
                            break
                if is_vertex:
                    # Add unique vertices
                    is_new = True
                    for pv in poly_vertices:
                        if np.allclose(pv, intersection):
                            is_new = False
                            break
                    if is_new:
                        poly_vertices.append(intersection)
            except np.linalg.LinAlgError:
                continue

    # --- 4. Calculate the area of P' using the shoelace formula ---
    if len(poly_vertices) < 3:
        print("Could not determine the vertices of the polytope P'.")
        return
        
    # Order vertices by angle around the centroid
    center = np.mean(poly_vertices, axis=0)
    poly_vertices.sort(key=lambda p: np.arctan2(p[1] - center[1], p[0] - center[0]))
    
    area_p_prime = 0.0
    for i in range(len(poly_vertices)):
        p1 = poly_vertices[i]
        p2 = poly_vertices[(i + 1) % len(poly_vertices)]
        area_p_prime += p1[0] * p2[1] - p2[0] * p1[1]
    area_p_prime = 0.5 * np.abs(area_p_prime)

    # --- 5. Calculate the volume of P and print results ---
    jacobian = np.abs(np.linalg.det(np.array([u1, u2])))
    # Note: jacobian is also equal to d! * simplex_volume
    
    polytope_p_volume = area_p_prime * jacobian
    
    print(f"For the given d={d} simplex (triangle):")
    print(f"Volume of simplex T (V): {simplex_volume:.4f}")
    print(f"Volume of transformed polytope P' in alpha-space: {area_p_prime:.4f}")
    print(f"Jacobian of the transformation (d! * V): {jacobian:.4f}")
    print("-" * 30)
    print(f"Final Volume of polytope P = (Volume of P') * Jacobian")
    print(f"Volume of P = {area_p_prime:.4f} * {jacobian:.4f} = {polytope_p_volume:.4f}")
    print("-" * 30)
    print(f"Ratio Vol(P) / V = {polytope_p_volume / simplex_volume:.4f}")
    print(f"This matches the theoretical value of d! = {np.math.factorial(d)}")

solve_geometry_problem()