import numpy as np
from scipy.spatial.transform import Rotation

def solve_icosahedron_problem():
    """
    Calculates and describes the shape of the water surface in a half-filled
    icosahedron tank standing on one of its faces.
    """
    # Step 1: Define Icosahedron Geometry
    # The golden ratio, phi
    PHI = (1 + np.sqrt(5)) / 2

    # Icosahedron vertices (these form three orthogonal golden rectangles)
    # The edge length for this set of vertices is 2.
    vertices = np.array([
        [-1,  PHI,  0], [ 1,  PHI,  0], [-1, -PHI,  0], [ 1, -PHI,  0],
        [ 0, -1,  PHI], [ 0,  1,  PHI], [ 0, -1, -PHI], [ 0,  1, -PHI],
        [PHI,  0, -1], [PHI,  0,  1], [-PHI,  0, -1], [-PHI,  0,  1]
    ])

    # Icosahedron faces are defined by indices into the vertices array.
    faces = np.array([
        [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
        [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
        [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
        [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1]
    ])

    # Define the 30 unique edges from the faces for intersection checks.
    edges = set()
    for face in faces:
        for i in range(3):
            v1_idx, v2_idx = face[i], face[(i + 1) % 3]
            # Store as a sorted tuple to ensure each edge is added only once.
            edge = tuple(sorted((v1_idx, v2_idx)))
            edges.add(edge)
    edges = np.array(list(edges))

    # Step 2: Orient the Icosahedron
    # Pick the first face to be the base and calculate its normal vector.
    base_face_indices = faces[0]
    v1, v2, v3 = vertices[base_face_indices]
    normal = np.cross(v2 - v1, v3 - v1)
    normal /= np.linalg.norm(normal)

    # We want the base to be flat on the ground, so its normal should point down.
    target_normal = np.array([0, 0, -1.0])

    # Calculate the rotation required to align the face's normal with the target normal.
    rotation, _ = Rotation.align_vectors([target_normal], [normal])
    rotated_vertices = rotation.apply(vertices)

    # Step 3: Translate to z=0
    # After rotation, all vertices of the base face have the same z-coordinate.
    # Translate the whole icosahedron so this base rests at z=0.
    base_z = rotated_vertices[base_face_indices[0], 2]
    translated_vertices = rotated_vertices - np.array([0, 0, base_z])

    # Step 4: Find the Half-Fill Level (Mid-height)
    # The total height is the difference between the max and min z-coordinates.
    z_min = np.min(translated_vertices[:, 2])
    z_max = np.max(translated_vertices[:, 2])
    height = z_max - z_min
    mid_height = height / 2

    # Step 5 & 6: Find intersection points with the mid-height plane.
    intersection_points = []
    for v1_idx, v2_idx in edges:
        p1 = translated_vertices[v1_idx]
        p2 = translated_vertices[v2_idx]

        # Check if the edge crosses the horizontal plane at mid_height.
        if (p1[2] < mid_height and p2[2] > mid_height) or \
           (p1[2] > mid_height and p2[2] < mid_height):
            # Calculate the intersection point using linear interpolation.
            t = (mid_height - p1[2]) / (p2[2] - p1[2])
            intersection_point = p1 + t * (p2 - p1)
            intersection_points.append(intersection_point)

    # Step 7: Analyze and Report the Resulting Polygon
    # Round points to handle floating point inaccuracies and find unique vertices.
    unique_points = np.unique(np.array(intersection_points).round(decimals=5), axis=0)
    num_vertices = len(unique_points)

    shape_description = f"a polygon with {num_vertices} vertices"
    if num_vertices == 6:
        shape_description = "a regular hexagon"

    print("An icosahedron is a 20-sided polyhedron where every face is an equilateral triangle.")
    print("When it stands on one face, its base is a triangle and its top is a parallel triangle.")
    print("Due to the icosahedron's symmetry, the half-volume level is exactly at its mid-height.")
    print("\nBy calculating the cross-section of the icosahedron at this mid-height plane:")
    print(f"The calculation found a polygon with {num_vertices} vertices.")
    print(f"\nConclusion: The shape of the water surface will be {shape_description}.")

    # Output properties of the resulting shape, per the user's request.
    if num_vertices == 6:
        # Sort points by angle to properly connect them and calculate side lengths.
        angles = np.arctan2(unique_points[:, 1], unique_points[:, 0])
        sorted_indices = np.argsort(angles)
        sorted_points = unique_points[sorted_indices]

        # Calculate the distance between consecutive vertices (side lengths).
        side_lengths = np.linalg.norm(np.roll(sorted_points, -1, axis=0) - sorted_points, axis=1)
        avg_side_length = np.mean(side_lengths)
        
        # Calculate the distance from the center to each vertex (radius).
        radii = np.linalg.norm(unique_points[:, :2], axis=1)
        avg_radius = np.mean(radii)

        print("\nProperties of the resulting hexagon:")
        print(f"  - The shape is centered at (x=0.0, y=0.0) on the plane z={mid_height:.4f}.")
        print(f"  - Its average side length is {avg_side_length:.4f} units.")
        print(f"  - Its average radius (distance from center to vertex) is {avg_radius:.4f} units.")
        if np.isclose(avg_side_length, avg_radius):
            print("  - As expected for a regular hexagon, the side length is equal to the radius.")
        
        print("\nThe shape is defined by its vertices. The coordinates of the hexagon's vertices are:")
        for i, p in enumerate(sorted_points):
            print(f"  Vertex {i+1}: (x={p[0]:.4f}, y={p[1]:.4f}, z={p[2]:.4f})")

if __name__ == '__main__':
    solve_icosahedron_problem()