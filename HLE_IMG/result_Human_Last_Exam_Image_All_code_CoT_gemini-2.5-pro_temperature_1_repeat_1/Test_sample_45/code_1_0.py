import math

def find_rotational_symmetry():
    """
    This function explains the logic for finding the rotational symmetry of the tiling
    and prints the result.
    """
    print("Step 1: Understanding Rotational Symmetry")
    print("A tiling has n-fold rotational symmetry if it remains unchanged after a rotation of 360/n degrees around a central point.")
    print("We need to find the highest possible integer 'n' > 1 for the given pattern.\n")

    print("Step 2: Analyzing the Tiling Pattern")
    print("Let's choose the center of one of the dark blue hexagons as our point of rotation.")
    print("By observing the pattern, we can see that each blue hexagon is surrounded by a ring of 6 yellow pentagons.\n")

    print("Step 3: Determining the Order of Symmetry")
    print("This structure of 6 identical shapes arranged symmetrically around a central point indicates that the pattern will look the same after being rotated 6 times within a full 360-degree turn.")
    order_n = 6
    print(f"The number of repeating units around the center is {order_n}. Therefore, the order of rotational symmetry is {order_n}.\n")

    print("Step 4: Calculating the Angle of Rotation")
    total_degrees = 360
    angle = total_degrees / order_n
    print(f"The smallest angle of rotation that leaves the pattern unchanged can be calculated.")
    print(f"The equation is: {total_degrees} / {order_n} = {int(angle)}")
    print(f"This confirms a {order_n}-fold rotational symmetry.\n")

    print("Final Answer:")
    print("The rotational symmetry of the tiling is:")
    print(order_n)

find_rotational_symmetry()