def find_rotational_symmetry():
    """
    This function analyzes the rotational symmetry of the provided tiling.
    """
    print("1. Rotational symmetry is determined by finding a point around which the tiling can be rotated and remain unchanged.")
    print("2. Let's choose the center of one of the prominent blue hexagons as our center of rotation.")
    print("3. By observing the pattern, we can see that there are 6 identical arrangements of yellow and orange tiles surrounding each blue hexagon.")
    
    # The order 'n' is the number of times the pattern repeats in a 360-degree rotation.
    n_fold_symmetry = 6
    
    print(f"4. This means the pattern repeats {n_fold_symmetry} times in a full circle.")
    print(f"5. Therefore, the tiling has a {n_fold_symmetry}-fold rotational symmetry.")
    
    # As requested, output the equation for the angle of rotation, including all numbers.
    total_degrees = 360
    rotation_angle = total_degrees / n_fold_symmetry
    
    print("\nThe equation for the angle of one rotational step is:")
    print(f"{total_degrees} / {n_fold_symmetry} = {int(rotation_angle)}")
    
    print(f"\nThe rotational symmetry of the tiling is {n_fold_symmetry}.")

find_rotational_symmetry()