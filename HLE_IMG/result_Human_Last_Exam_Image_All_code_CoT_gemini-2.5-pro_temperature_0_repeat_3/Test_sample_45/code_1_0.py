def find_rotational_symmetry():
    """
    This function analyzes the rotational symmetry of the provided tiling.
    """
    print("Analyzing the rotational symmetry of the tiling...")

    # Step 1: Identify potential centers of rotation and their orders.
    # By inspecting the pattern, we can identify centers of rotation at the center of the hexagons.
    order_at_hexagon = 6
    angle_at_hexagon = 360 / order_at_hexagon

    print(f"A rotation around the center of a hexagon by {int(angle_at_hexagon)} degrees leaves the pattern unchanged.")
    print(f"This corresponds to a {order_at_hexagon}-fold rotational symmetry.")
    print(f"Calculation: 360 / {order_at_hexagon} = {int(angle_at_hexagon)}")

    # We can also observe apparent symmetry at the center of the squares.
    order_at_square = 4
    angle_at_square = 360 / order_at_square
    print(f"\nA rotation around the center of a square by {int(angle_at_square)} degrees appears to leave the pattern unchanged.")
    print(f"This corresponds to a {order_at_square}-fold rotational symmetry.")
    print(f"Calculation: 360 / {order_at_square} = {int(angle_at_square)}")

    # Step 2: Determine the overall rotational symmetry of the tiling.
    # The rotational symmetry of the entire pattern is defined by the highest order of symmetry found.
    highest_order = max(order_at_hexagon, order_at_square)

    print("\nThe overall rotational symmetry of a pattern is the highest order of rotational symmetry it possesses.")
    print(f"The highest order found is {highest_order}.")

find_rotational_symmetry()