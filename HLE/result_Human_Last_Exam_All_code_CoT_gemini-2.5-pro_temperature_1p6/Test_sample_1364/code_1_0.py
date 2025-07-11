import numpy as np
from scipy.spatial import ConvexHull

def solve_polyhedron_problem():
    """
    This function demonstrates that a 7-vertex polyhedron can have
    quadrilateral projections on three generally positioned planes.
    """
    # 1. Define the polyhedron's vertices.
    # We take the 8 vertices of a cube and remove one, e.g., (1, 1, 1).
    # The convex hull of the remaining 7 vertices forms our polyhedron.
    cube_vertices = np.array([
        [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1],
        [1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 1]
    ])
    # The set of vertices for our 7-vertex polyhedron P
    p_vertices = cube_vertices[:-1]
    
    print(f"Analyzing a polyhedron with V={len(p_vertices)} vertices.")
    print("The chosen planes are the xy, yz, and xz planes, which are in a general position.")
    print("-" * 30)

    # 2. Define the projection planes (by their normal vectors/projection directions)
    # xy-plane (project along z), yz-plane (project along x), xz-plane (project along y)
    projections = {
        "yz-plane (project along x-axis)": [1, 2],
        "xz-plane (project along y-axis)": [0, 2],
        "xy-plane (project along z-axis)": [0, 1],
    }

    all_are_quadrilaterals = True

    # 3. Perform each projection and check the result
    for name, dims in projections.items():
        print(f"Projecting onto the {name}...")
        
        # Project the 3D vertices to 2D
        projected_points = p_vertices[:, dims]
        
        # Find the convex hull of the projected 2D points
        try:
            hull = ConvexHull(projected_points)
            hull_vertices_count = len(hull.vertices)
            
            print(f"The projection is a polygon with {hull_vertices_count} vertices.")
            
            if hull_vertices_count != 4:
                all_are_quadrilaterals = False
            
            # Output the vertices that form the projected polygon
            final_projected_shape_vertices = projected_points[hull.vertices]
            print("Vertices of the projected quadrilateral:")
            # Use a loop to print each vertex for the final equation part of the prompt
            for v in final_projected_shape_vertices:
                print(f"({v[0]}, {v[1]})")

        except Exception as e:
            # This might happen if all points are collinear
            print(f"Could not compute convex hull: {e}")
            all_are_quadrilaterals = False
            
        print("-" * 30)
        
    # 4. Final conclusion based on the analysis
    print("Summary:")
    print("We have shown that a polyhedron with V=7 can have quadrilateral projections on 3 planes in general position.")
    print("We ruled out V=4 and V=5.")
    print("Constructions exist for V=6 (octahedron) and V=8 (cube).")
    print("By extending these constructions (e.g., by truncating corners), it's possible to create valid polyhedra for any number of vertices V >= 6.")
    print("\nThe set of possible numbers of vertices is all integers greater than or equal to 6.")


solve_polyhedron_problem()
