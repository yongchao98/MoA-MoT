import math

def solve_symmetry():
    """
    This function determines the rotational symmetry of the tiling shown in the image.
    """
    # A full circle rotation is 360 degrees.
    full_rotation_degrees = 360

    # By observing the pattern, we can identify centers of rotation.
    # Let's focus on the center of one of the dark blue hexagonal shapes.
    # The pattern of surrounding yellow squares and orange rhombuses is repeated
    # multiple times around this center.
    # Counting these repetitions reveals how many times the pattern can be rotated
    # to match itself in a full 360-degree turn.
    # Upon inspection, we can see 6 identical arrangements around the center.
    number_of_repetitions = 6

    # The order of rotational symmetry is equal to this number of repetitions.
    order_of_symmetry = number_of_repetitions

    # The minimum angle of rotation that leaves the tiling unchanged can be calculated
    # by dividing the full circle by the order of symmetry.
    minimum_rotation_angle = full_rotation_degrees / order_of_symmetry

    print("Analyzing the rotational symmetry of the tiling:")
    print(f"1. We identify a center of rotation at the middle of a blue hexagonal shape.")
    print(f"2. We count the number of times the surrounding pattern repeats itself around this center. The count is {number_of_repetitions}.")
    print(f"3. The order of rotational symmetry is equal to this number of repetitions.")
    print(f"4. The minimum rotation angle is {full_rotation_degrees} degrees / {number_of_repetitions} = {int(minimum_rotation_angle)} degrees.")
    print(f"5. The order of symmetry is the number of times the pattern maps onto itself in a full 360-degree turn.")
    print(f"The calculation for the order is: {full_rotation_degrees} / {int(minimum_rotation_angle)} = {order_of_symmetry}")
    print(f"Therefore, the tiling has a {order_of_symmetry}-fold rotational symmetry.")


solve_symmetry()