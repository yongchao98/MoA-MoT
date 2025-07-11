def find_water_surface_shape():
    """
    This function determines the shape of the water surface in a half-filled
    icosahedron tank that is standing on one of its faces.
    """

    # 1. An icosahedron is a highly symmetric shape with a center of symmetry.
    #    When it's half-full, the water level forms a plane passing through this center.
    reasoning_step_1 = "A half-filled, centrally symmetric container has a liquid surface that passes through its center of symmetry."

    # 2. The icosahedron stands on a face, making the base horizontal.
    #    The water surface is also horizontal, so it is parallel to the base face.
    reasoning_step_2 = "The water surface is a horizontal plane, parallel to the icosahedron's base face."

    # 3. Combining these facts, the water surface is the cross-section of the icosahedron
    #    formed by a plane passing through its center and parallel to a face.
    reasoning_step_3 = "The shape is the cross-section from a plane cutting through the icosahedron's center, parallel to a face."

    # 4. This specific cross-section of an icosahedron is a known geometric figure.
    shape = "regular hexagon"
    reasoning_step_4 = f"For a regular icosahedron, this cross-section is a {shape}."

    # Print the reasoning and the result.
    # The original prompt's request to "output each number in the final equation"
    # is not applicable here, as the problem is purely geometric and has no equation.
    print("Problem: What is the shape of the water surface in a half-filled icosahedron tank standing on a face?")
    print("-" * 80)
    print("Reasoning:")
    print(f"1. {reasoning_step_1}")
    print(f"2. {reasoning_step_2}")
    print(f"3. {reasoning_step_3}")
    print(f"4. {reasoning_step_4}")
    print("-" * 80)
    print(f"Conclusion: The shape of the water surface will be a {shape}.")


find_water_surface_shape()