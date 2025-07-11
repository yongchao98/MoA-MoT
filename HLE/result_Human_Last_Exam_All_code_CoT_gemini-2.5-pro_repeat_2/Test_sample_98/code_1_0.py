def get_water_surface_shape():
    """
    Analyzes the geometry of a half-filled icosahedron tank to determine
    the shape of the water's surface.
    """

    # Step 1: Define the physical setup.
    # The tank is an icosahedron, which has 20 identical equilateral triangle faces.
    # It stands on one face, meaning that face is the horizontal base.
    # Due to the icosahedron's symmetry, there is a parallel top face.
    print("Step 1: The icosahedron tank stands on a triangular face, which is its base.")

    # Step 2: Understand the "half-filled" condition.
    # When the tank is half-filled, the water occupies exactly 50% of the volume.
    # The surface of the water is a horizontal plane.
    # For the volume to be precisely half, this plane must pass through the center
    # of symmetry of the icosahedron.
    print("Step 2: Being half-full means the water surface is a horizontal plane passing through the center of the icosahedron.")

    # Step 3: Identify which faces the water surface intersects.
    # This central plane is parallel to the top and bottom faces.
    # It is located at the "equator" of the icosahedron in this orientation.
    # This equatorial region, or "belt," consists of 6 triangular faces.
    # The water surface plane cuts through each of these 6 faces.
    print("Step 3: This central plane intersects the 6 triangular faces forming the 'belt' of the icosahedron.")

    # Step 4: Determine the shape from the intersections.
    # The intersection of a plane with 6 faces forms a polygon with 6 sides.
    # Therefore, the shape of the water surface is a hexagon.
    print("Step 4: The intersection with 6 faces results in a 6-sided polygon (a hexagon).")

    # Step 5: Use symmetry to determine the type of hexagon.
    # The icosahedron has a 3-fold rotational symmetry axis passing vertically
    # through the centers of the base and top faces. The horizontal water surface
    # is perpendicular to this axis, so the resulting hexagon must also have
    # 3-fold rotational symmetry.
    # The plane also passes through the icosahedron's center of inversion, which
    # means the hexagon must be centrally symmetric.
    # A hexagon with both 3-fold rotational symmetry and central symmetry is
    # required to be a regular hexagon.
    print("Step 5: The symmetries of the icosahedron require this hexagon to be regular.")

    # Final Conclusion
    final_shape = "A regular hexagon"
    print("\nConclusion: The shape of the water surface is a regular hexagon.")


if __name__ == "__main__":
    get_water_surface_shape()