def calculate_euler_characteristic_2d(vertices, edges, faces):
    """
    Calculates the Euler characteristic for a 2D simplicial complex (triangulated surface).

    The Euler characteristic is defined as V - E + F, where:
    - V is the number of vertices (points)
    - E is the number of edges (lines connecting vertices)
    - F is the number of faces (triangles)

    Args:
        vertices (list or set): A collection of vertices.
        edges (list or set): A collection of edges.
        faces (list or set): A collection of faces.

    Returns:
        int: The Euler characteristic of the complex.
    """
    num_vertices = len(vertices)
    num_edges = len(edges)
    num_faces = len(faces)
    
    chi = num_vertices - num_edges + num_faces
    
    print(f"Number of Vertices (V): {num_vertices}")
    print(f"Number of Edges (E): {num_edges}")
    print(f"Number of Faces (F): {num_faces}")
    print(f"Euler Characteristic (V - E + F) = {num_vertices} - {num_edges} + {num_faces} = {chi}")
    
    return chi

# Example: A triangulated annulus (a ring shape).
# This surface is non-compact (if viewed as open) or has a boundary.
# Its Euler characteristic is 0.
# We can construct it from 6 vertices and 6 triangular faces.
# Vertices: 3 on an inner ring, 3 on an outer ring.
# Faces: connecting inner and outer vertices.

# Define the components of a triangulated annulus
annulus_vertices = {'v_in1', 'v_in2', 'v_in3', 'v_out1', 'v_out2', 'v_out3'}
annulus_edges = {
    # Inner edges
    ('v_in1', 'v_in2'), ('v_in2', 'v_in3'), ('v_in3', 'v_in1'),
    # Outer edges
    ('v_out1', 'v_out2'), ('v_out2', 'v_out3'), ('v_out3', 'v_out1'),
    # Radial edges
    ('v_in1', 'v_out1'), ('v_in2', 'v_out2'), ('v_in3', 'v_out3'),
    ('v_in1', 'v_out2'), ('v_in2', 'v_out3'), ('v_in3', 'v_out1') # These are the diagonals of the quadrilaterals
}
annulus_faces = {
    ('v_in1', 'v_out1', 'v_out2'), ('v_in1', 'v_in2', 'v_out2'),
    ('v_in2', 'v_out2', 'v_out3'), ('v_in2', 'v_in3', 'v_out3'),
    ('v_in3', 'v_out3', 'v_out1'), ('v_in3', 'v_in1', 'v_out1')
}

print("Calculating Euler characteristic for a triangulated annulus:")
calculate_euler_characteristic_2d(annulus_vertices, annulus_edges, annulus_faces)

# The result being 0 is consistent with the fact that such a surface admits
# a non-vanishing vector field, which is related to the existence of a section
# for the configuration space fibration.