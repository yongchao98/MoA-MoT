def solve_rotational_symmetry():
    """
    This function calculates and explains the rotational symmetry of the provided tiling.
    """
    
    # A full rotation is 360 degrees.
    total_degrees = 360

    # By observing the tiling, we can identify a center of rotation.
    # Let's choose the center of a dark blue hexagon.
    # We then count the number of times the surrounding pattern repeats itself in a full circle.
    # Around the central hexagon, there are 6 identical arrangements of squares and rhombuses.
    order_of_symmetry = 6

    # The angle of rotation is the total degrees divided by the order of symmetry.
    angle_of_rotation = total_degrees / order_of_symmetry

    print("Analyzing the rotational symmetry of the tiling:")
    print(f"1. A full circle is {total_degrees} degrees.")
    print(f"2. We choose a point of rotation, like the center of a hexagon.")
    print(f"3. We count the number of identical repeating patterns around this point. There are {order_of_symmetry} repeats.")
    print(f"4. The order of rotational symmetry is therefore {order_of_symmetry}.")
    print(f"5. The smallest angle of rotation that leaves the pattern unchanged is calculated as:")
    print(f"   {total_degrees} / {order_of_symmetry} = {int(angle_of_rotation)} degrees.")
    print("\nConclusion: The rotational symmetry of the tiling is " + str(order_of_symmetry) + ".")

solve_rotational_symmetry()
