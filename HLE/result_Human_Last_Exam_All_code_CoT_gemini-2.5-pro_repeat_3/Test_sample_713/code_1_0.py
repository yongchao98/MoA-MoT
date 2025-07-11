import numpy as np

def get_polytope_area_for_triangle(vertices):
    """
    Calculates the area of the polytope P for a given triangle T in 2D.
    
    Args:
        vertices (list of tuples): A list of 3 tuples, where each tuple is a 2D vertex of the triangle.
        
    Returns:
        tuple: A tuple containing the area of the polytope P, the area of the triangle T, and their ratio.
    """
    v = [np.array(p) for p in vertices]
    v0, v1, v2 = v[0], v[1], v[2]

    # 1. Calculate the area of the triangle T
    triangle_area = 0.5 * np.abs(np.linalg.det(np.array([v1 - v0, v2 - v0])))

    # 2. Define the hyperplanes (lines) for the polytope P
    edges = [(v0, v1), (v0, v2), (v1, v2)]
    slabs = []
    lines = []

    for p1, p2 in edges:
        normal = p2 - p1
        # Ensure normal is not a zero vector
        if np.linalg.norm(normal) == 0:
            continue
        
        c1 = np.dot(normal, p1)
        c2 = np.dot(normal, p2)
        
        # A slab is defined by a normal vector and two constants c_min, c_max
        # such that c_min <= normal.x <= c_max
        slabs.append({'normal': normal, 'c_min': min(c1, c2), 'c_max': max(c1, c2)})
        
        # A line is defined by ax+by=c
        lines.append({'normal': normal, 'c': c1})
        lines.append({'normal': normal, 'c': c2})

    # 3. Find the vertices of the polytope P
    poly_vertices = []
    # Find intersections of all pairs of lines
    for i in range(len(lines)):
        for j in range(i + 1, len(lines)):
            line1 = lines[i]
            line2 = lines[j]
            
            A = np.array([line1['normal'], line2['normal']])
            b = np.array([line1['c'], line2['c']])
            
            # If lines are parallel, they don't have a unique intersection
            if np.abs(np.linalg.det(A)) < 1e-9:
                continue
                
            # Solve for the intersection point
            try:
                intersection_point = np.linalg.solve(A, b)
            except np.linalg.LinAlgError:
                continue

            # Check if the point is within all defined slabs
            is_vertex = True
            for slab in slabs:
                val = np.dot(slab['normal'], intersection_point)
                if not (slab['c_min'] - 1e-9 <= val <= slab['c_max'] + 1e-9):
                    is_vertex = False
                    break
            
            if is_vertex:
                # Add vertex if it's not already in the list (within tolerance)
                is_new = True
                for pv in poly_vertices:
                    if np.linalg.norm(pv - intersection_point) < 1e-9:
                        is_new = False
                        break
                if is_new:
                    poly_vertices.append(intersection_point)

    # 4. Calculate the area of the convex polygon P using the shoelace formula
    if len(poly_vertices) < 3:
        return 0, triangle_area, 0

    # Sort vertices by angle around the centroid
    center = np.mean(poly_vertices, axis=0)
    poly_vertices.sort(key=lambda p: np.arctan2(p[1] - center[1], p[0] - center[0]))
    
    x = np.array([p[0] for p in poly_vertices])
    y = np.array([p[1] for p in poly_vertices])
    
    polytope_area = 0.5 * np.abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))

    # 5. Compute the ratio
    ratio = polytope_area / triangle_area if triangle_area > 0 else 0
    
    return polytope_area, triangle_area, ratio

def main():
    """
    Main function to run the demonstration for d=2.
    """
    print("Demonstrating the result for d=2.")
    print("The theoretical ratio Vol(P)/Vol(T) should be 2! = 2.")
    print("-" * 50)

    # Case 1: A right-angled triangle
    right_triangle_vertices = [(0, 0), (3, 0), (0, 4)]
    p_area, t_area, ratio = get_polytope_area_for_triangle(right_triangle_vertices)
    
    print("Case 1: Right-angled triangle with vertices at (0,0), (3,0), (0,4)")
    print(f"Area of Triangle (T): {t_area:.4f}")
    print(f"Area of Polytope (P): {p_area:.4f}")
    print(f"Ratio Area(P) / Area(T) = {ratio:.4f}")
    print("Final Equation: Area(P) / Area(T) = 2")
    print("-" * 50)
    
    # Case 2: An equilateral triangle
    s = 4
    h = s * np.sqrt(3) / 2
    equilateral_triangle_vertices = [(0, 0), (s, 0), (s/2, h)]
    p_area, t_area, ratio = get_polytope_area_for_triangle(equilateral_triangle_vertices)
    
    print(f"Case 2: Equilateral triangle with side length {s}")
    print(f"Area of Triangle (T): {t_area:.4f}")
    print(f"Area of Polytope (P): {p_area:.4f}")
    print(f"Ratio Area(P) / Area(T) = {ratio:.4f}")
    print("Final Equation: Area(P) / Area(T) = 2")
    print("-" * 50)
    
    print("The results for d=2 are consistent with the general formula.")
    print("\nFor a non-degenerate simplex T in d-dimensional space with volume V,")
    print("the volume of the corresponding polytope P is given by the equation:")
    print("Vol(P) = d! * V")
    
if __name__ == "__main__":
    main()
