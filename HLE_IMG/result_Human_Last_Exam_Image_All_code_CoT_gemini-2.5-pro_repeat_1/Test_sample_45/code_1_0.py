def find_rotational_symmetry():
    """
    Analyzes the rotational symmetry of the tiling shown in the image.

    The tiling is composed of hexagons, squares, and rhombuses. We look for the
    point of highest rotational symmetry.

    1. The center of each light yellow square is a center of 4-fold rotation.
       A rotation of 360 / 4 = 90 degrees leaves the pattern unchanged.

    2. The center of each dark blue hexagon is a center of 6-fold rotation.
       A rotation of 360 / 6 = 60 degrees leaves the pattern unchanged.

    The overall rotational symmetry of a pattern is determined by the highest
    order of rotation found. In this case, it is 6.
    """
    highest_order_of_symmetry = 6
    angle_of_rotation = 360 / highest_order_of_symmetry

    print(f"The analysis shows the presence of different rotation centers.")
    print(f"The highest order of rotation is found at the center of the dark blue hexagons.")
    print(f"The tiling remains unchanged after a rotation of 360 / {highest_order_of_symmetry} = {int(angle_of_rotation)} degrees.")
    print(f"Therefore, the rotational symmetry of the tiling is {highest_order_of_symmetry}.")

find_rotational_symmetry()