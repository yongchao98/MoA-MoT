import math

def generate_polyhedron_vertices(F):
    """
    Generates and prints the vertices of a convex polyhedron with F vertices
    that satisfies the problem's condition.

    Args:
        F (int): The desired number of vertices, must be >= 4.
    """
    if not isinstance(F, int) or F < 4:
        print("The number of vertices F must be an integer >= 4.")
        return

    print(f"--- Generating vertices for a polyhedron with F = {F} vertices ---")
    vertices = []
    if F == 4:
        # For F=4, we use a specific tetrahedron whose projection on xy, yz, xz planes is a square.
        print("This is a tetrahedron.")
        vertices.extend([
            (1, 1, 1),
            (1, -1, -1),
            (-1, 1, -1),
            (-1, -1, 1)
        ])
    else:
        # For F >= 5, we use the general construction: a square base with a cap of k = F-4 vertices.
        # The number of vertices in the cap is k = F - 4.
        k = F - 4
        print(f"This polyhedron is constructed from a square base and a cap of {k} vertex/vertices.")
        
        # 1. Add the 4 vertices of the square base S
        base_vertices = [
            (1.0, 1.0, 0.0),
            (1.0, -1.0, 0.0),
            (-1.0, -1.0, 0.0),
            (-1.0, 1.0, 0.0)
        ]
        vertices.extend(base_vertices)

        # 2. Add the k vertices of the cap U
        if k > 0:
            radius = 0.5  # Radius for placing cap vertices, must be < 1
            for j in range(k):
                angle = 2 * math.pi * j / k
                x = radius * math.cos(angle)
                y = radius * math.sin(angle)
                z = 1.0
                vertices.append((round(x, 5), round(y, 5), round(z, 5)))

    # Print the final list of vertices
    print("Vertices:")
    for i, v in enumerate(vertices):
        print(f"  v{i+1}: {v}")
    print("-" * (35 + len(str(F))))


# --- Main execution ---
# Generate examples for a few possible numbers of vertices.
generate_polyhedron_vertices(4)
generate_polyhedron_vertices(5)
generate_polyhedron_vertices(7)
generate_polyhedron_vertices(8)
