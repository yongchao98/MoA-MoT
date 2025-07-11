def solve_rotational_symmetry():
    """
    This function determines the rotational symmetry of the tiling shown in the image.
    """
    # The rotational symmetry is determined by the highest order of rotation found in the pattern.
    # We identify that the centers of the dark blue hexagons provide the highest order.
    # A hexagon has 6 sides, which corresponds to a 6-fold symmetry.
    
    num_sides_of_central_polygon = 6
    full_rotation = 360
    
    # The angle of rotation for one step of the symmetry operation.
    rotation_angle = full_rotation / num_sides_of_central_polygon
    
    print("The rotational symmetry of a tiling is the highest order of rotation that leaves the pattern unchanged.")
    print("In this tiling, the highest order of symmetry is centered on the hexagons.")
    print(f"A hexagon has {num_sides_of_central_polygon} sides.")
    print("The calculation for the rotation angle is:")
    print(f"{full_rotation} (degrees in a full circle) / {num_sides_of_central_polygon} (sides of hexagon) = {int(rotation_angle)} degrees.")
    print(f"Since the pattern repeats {num_sides_of_central_polygon} times in a 360-degree rotation, the rotational symmetry is {num_sides_of_central_polygon}-fold.")
    
    # The final answer is the order of the symmetry.
    final_answer = num_sides_of_central_polygon
    print(f"\nThe rotational symmetry is: {final_answer}")

solve_rotational_symmetry()