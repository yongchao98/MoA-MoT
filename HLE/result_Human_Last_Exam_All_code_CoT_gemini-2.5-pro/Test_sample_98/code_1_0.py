def solve_water_surface_shape():
    """
    This script determines the shape of the water surface in a half-filled
    icosahedron tank standing on one of its faces by analyzing its geometry.
    """

    # 1. Analyze the structure of an icosahedron standing on a face.
    # Its 12 vertices are arranged in four horizontal layers of three vertices each.
    vertices_in_bottom_layer = 3
    vertices_in_lower_middle_layer = 3
    vertices_in_upper_middle_layer = 3
    vertices_in_top_layer = 3
    total_vertices = sum([
        vertices_in_bottom_layer,
        vertices_in_lower_middle_layer,
        vertices_in_upper_middle_layer,
        vertices_in_top_layer
    ])

    # 2. The water surface in a half-filled tank is a horizontal plane
    # passing through the center of the icosahedron. This plane lies
    # exactly between the two middle layers of vertices.

    # 3. This central plane cuts the edges that connect the lower-middle
    # layer of vertices to the upper-middle layer. In this orientation,
    # there are exactly 6 such edges forming the "waist" of the icosahedron.
    number_of_edges_cut = 6

    # 4. The number of edges cut by the plane equals the number of vertices
    # (and therefore sides) of the cross-sectional shape.
    number_of_sides = number_of_edges_cut

    # 5. Due to the icosahedron's rotational and reflection symmetries,
    # all sides of the resulting polygon are of equal length and all
    # interior angles are equal. A six-sided polygon with these properties
    # is a regular hexagon.
    shape_name = "regular hexagon"

    print("Problem: What is the shape of the water surface in a half-filled icosahedron tank standing on a face?")
    print("\nReasoning:")
    print(f"1. An icosahedron has {total_vertices} vertices. When stood on a face, these vertices form 4 layers of {vertices_in_bottom_layer} each.")
    print("2. A half-full tank means the water surface is a central plane slicing the icosahedron in half.")
    print(f"3. This plane intersects the {number_of_edges_cut} edges that connect the two middle layers of vertices.")
    print(f"4. The resulting cross-section is a polygon with {number_of_sides} sides.")
    print(f"5. Due to the perfect symmetry of the icosahedron, this shape is a {shape_name}.")

    print("\nFinal Answer:")
    print(f"The shape of the water surface will be a {shape_name}.")

solve_water_surface_shape()