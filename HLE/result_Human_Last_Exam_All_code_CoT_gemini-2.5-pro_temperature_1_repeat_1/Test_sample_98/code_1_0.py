def solve_tank_shape():
    """
    This function explains and prints the shape of the water surface
    in a half-filled icosahedron tank standing on a face.
    """

    # The shape is determined by the cross-section of the icosahedron
    # at the mid-plane between the top and bottom faces.
    # Based on the geometry and symmetry of the icosahedron, this
    # cross-section is a well-known shape.

    # Number of sides of the resulting polygon
    # This corresponds to the number of edges the mid-plane intersects.
    number_of_sides = 6

    # Description of the shape
    # Due to the high degree of symmetry, the hexagon is regular.
    shape_description = "regular hexagon"

    print("The shape of the water surface is a polygon with {} sides.".format(number_of_sides))
    print("Specifically, due to the symmetries of the icosahedron, the shape is a {}.".format(shape_description))

solve_tank_shape()