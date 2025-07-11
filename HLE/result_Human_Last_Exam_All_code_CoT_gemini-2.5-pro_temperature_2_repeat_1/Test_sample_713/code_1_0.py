import numpy as np

def triangle_polytope_area(v0, v1, v2):
    """
    Calculates the area of a triangle and its associated polytope P for d=2.
    
    Args:
      v0, v1, v2: Numpy arrays representing the vertices of the triangle.
      
    Returns:
      A tuple containing (Area(T), Area(P), Ratio Area(P)/Area(T)).
    """
    # 1. Calculate Triangle Area (V)
    # The area of a triangle with vertices (v0, v1, v2) can be calculated
    # using the determinant of the matrix of edge vectors.
    mat = np.array([v1 - v0, v2 - v0])
    V = 0.5 * np.abs(np.linalg.det(mat))

    # 2. Define the lines bounding the polytope P
    # P is bounded by 6 lines. Each pair of lines corresponds to an edge.
    edges = [(v0, v1), (v1, v2), (v2, v0)]
    lines = []
    normals = []
    for p1, p2 in edges:
        edge_vec = p2 - p1
        
        # Hyperplane through p1
        # (x - p1) . edge_vec = 0 => x . edge_vec = p1 . edge_vec
        c1 = np.dot(p1, edge_vec)
        lines.append((edge_vec, c1))
        
        # Hyperplane through p2
        # (x - p2) . edge_vec = 0 => x . edge_vec = p2 . edge_vec
        c2 = np.dot(p2, edge_vec)
        lines.append((edge_vec, c2))
        
        # Store normals for sorting
        normals.append(edge_vec)
        normals.append(-edge_vec)
        
    # Sort the lines by the angle of their normal vectors to find adjacent lines
    # This helps in finding the vertices of the convex hexagon in order.
    angles = [np.arctan2(n[1], n[0]) for n in [l[0] for l in lines]]
    sorted_indices = np.argsort(angles)
    sorted_lines = [lines[i] for i in sorted_indices]

    # 3. Find vertices of the hexagon P
    # Vertices are intersections of adjacent lines
    hex_vertices = []
    for i in range(len(sorted_lines)):
        n1, c1 = sorted_lines[i]
        n2, c2 = sorted_lines[(i + 1) % len(sorted_lines)]
        
        # Solve the system of linear equations:
        # n1[0]*x + n1[1]*y = c1
        # n2[0]*x + n2[1]*y = c2
        A = np.array([n1, n2])
        b = np.array([c1, c2])
        
        # Check if normals are parallel. If so, this pair doesn't form a vertex.
        if np.abs(np.linalg.det(A)) < 1e-9:
            continue
            
        vertex = np.linalg.solve(A, b)
        
        # To avoid duplicate vertices due to floating point inaccuracies
        is_duplicate = False
        for hv in hex_vertices:
            if np.linalg.norm(vertex - hv) < 1e-9:
                is_duplicate = True
                break
        if not is_duplicate:
            hex_vertices.append(vertex)

    # 4. Calculate hexagon area using Shoelace formula
    x_coords = [v[0] for v in hex_vertices]
    y_coords = [v[1] for v in hex_vertices]
    
    shoelace_sum1 = np.sum([x_coords[i] * y_coords[(i + 1) % len(hex_vertices)] for i in range(len(hex_vertices))])
    shoelace_sum2 = np.sum([y_coords[i] * x_coords[(i + 1) % len(hex_vertices)] for i in range(len(hex_vertices))])
    
    P_area = 0.5 * np.abs(shoelace_sum1 - shoelace_sum2)

    return V, P_area

def run_analysis(name, triangle_vertices):
    v0, v1, v2 = [np.array(v) for v in triangle_vertices]
    V, P_area = triangle_polytope_area(v0, v1, v2)
    print(f"--- Analysis for {name} Triangle ---")
    print(f"Vertices: v0={v0}, v1={v1}, v2={v2}")
    print(f"Simplex Volume V (Area): {V:.4f}")
    print(f"Polytope P Volume (Area): {P_area:.4f}")
    if V > 1e-9:
      print(f"Ratio Vol(P) / V: {P_area/V:.4f}")
    d = 2
    print(f"The theoretical minimum ratio for d={d} is d! = {np.math.factorial(d)}.")
    print("")

if __name__ == '__main__':
    # Case 1: Right-angled triangle (orthocentric)
    right_triangle = [(0, 0), (3, 0), (0, 4)]
    run_analysis("Right-angled", right_triangle)

    # Case 2: Equilateral triangle (orthocentric)
    equilateral_triangle = [(0, 0), (2, 0), (1, np.sqrt(3))]
    run_analysis("Equilateral", equilateral_triangle)
    
    # Case 3: Obtuse triangle (not orthocentric)
    obtuse_triangle = [(0, 0), (3, 0), (-1, 1)]
    run_analysis("Obtuse", obtuse_triangle)
    
    # Demonstrate the ratio can be very large
    very_obtuse_triangle = [(0,0), (3,0), (-2.9, 0.1)]
    run_analysis("Very Obtuse", very_obtuse_triangle)