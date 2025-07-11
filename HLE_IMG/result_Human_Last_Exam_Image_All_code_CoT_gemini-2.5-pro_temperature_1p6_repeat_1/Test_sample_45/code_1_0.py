def solve_rotational_symmetry():
    """
    This function explains and calculates the rotational symmetry of the provided tiling.
    """
    print("To find the rotational symmetry, we identify a point around which the tiling can be rotated and remain unchanged.")
    print("Let's choose the center of one of the dark blue hexagons as our point of rotation.")
    print("\nIf we rotate the entire pattern around this point, it will perfectly overlap with its original position 6 times in a full 360-degree turn.")
    print("This means the tiling has 6-fold rotational symmetry.")
    print("\nThe smallest positive angle of rotation for which the tiling looks the same is calculated by dividing 360 degrees by the order of symmetry.")

    total_degrees = 360
    order_of_symmetry = 6
    angle_of_rotation = total_degrees / order_of_symmetry

    print("\nCalculation:")
    print(f"{total_degrees} degrees / {order_of_symmetry} = {int(angle_of_rotation)} degrees")

solve_rotational_symmetry()