def solve_rotational_symmetry():
    """
    Calculates the rotational symmetry of the provided tiling.
    """
    total_degrees = 360
    # By visual inspection, the smallest angle of rotation that leaves the
    # tiling unchanged is 90 degrees. This rotation can be centered on
    # any of the yellow squares or the point where four orange rhombi meet.
    min_rotation_angle = 90

    # The order of rotational symmetry is the total degrees in a circle
    # divided by the minimum angle of rotation.
    order_of_symmetry = total_degrees / min_rotation_angle

    print(f"To find the rotational symmetry, we divide the total degrees in a circle (360) by the smallest angle of rotation that leaves the pattern unchanged.")
    print(f"From the image, we can see the smallest angle of rotation is {min_rotation_angle} degrees.")
    print(f"The calculation is:")
    # The prompt requests that each number in the final equation is printed.
    print(f"{total_degrees} / {min_rotation_angle} = {int(order_of_symmetry)}")
    print(f"Therefore, the tiling has a {int(order_of_symmetry)}-fold rotational symmetry.")

solve_rotational_symmetry()