def solve_icosahedron_problem():
    """
    Determines and prints the shape of the water surface in a half-filled icosahedron tank.
    """

    # An icosahedron is a regular polyhedron with 20 identical equilateral triangular faces.
    # It is standing on one face, which is horizontal.
    # By symmetry, there is a parallel face at the top.

    # When the tank is half-filled, the water level is at a horizontal plane
    # that bisects the volume. Due to symmetry, this plane is exactly halfway
    # between the top and bottom faces.

    # An icosahedron has 12 vertices.
    # - 3 vertices form the bottom face.
    # - 3 vertices form the top face.
    # - The remaining 6 vertices lie on a horizontal plane in the middle.

    # The water surface is the intersection of the mid-plane with the icosahedron.
    # This intersection is a polygon whose vertices are the 6 middle vertices.
    # This polygon is a hexagon.

    # The edges of this hexagon are also edges of the icosahedron.
    # Since all edges of a regular icosahedron are equal, the hexagon is equilateral.

    # Due to the 3-fold rotational symmetry of the icosahedron about the vertical axis,
    # the hexagon is also equiangular.

    # A hexagon that is both equilateral and equiangular is a regular hexagon.
    final_shape = "a regular hexagon"

    print("The problem asks for the shape of the water surface in a half-filled icosahedron tank standing on one of its faces.")
    print(f"Based on the symmetry and geometry of the icosahedron, the shape of the surface is {final_shape}.")

solve_icosahedron_problem()