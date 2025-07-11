def find_rotational_symmetry():
    """
    This function explains the reasoning for determining the rotational symmetry of the given tiling.
    """
    
    # The number of identical segments (the 'arms' of the star) around a central point.
    # By observing the star-like motifs in the tiling, we can count 5 identical arms.
    n_fold_symmetry = 5
    
    # The total degrees in a full circle.
    full_circle_degrees = 360
    
    # The smallest angle of rotation that leaves the pattern unchanged.
    angle_of_rotation = full_circle_degrees / n_fold_symmetry
    
    print("Step 1: Identify a central point for rotation in the tiling. The center of the star-like motifs is a suitable point.")
    print(f"Step 2: Count the number of times the pattern repeats in a full 360-degree rotation around this point. We can see {n_fold_symmetry} identical 'arms' forming the star.")
    print(f"Step 3: The order of rotational symmetry is equal to this count, which is {n_fold_symmetry}.")
    print(f"Step 4: This means a rotation of {full_circle_degrees} / {n_fold_symmetry} = {int(angle_of_rotation)} degrees leaves the tiling unchanged.")
    print("\nTherefore, the rotational symmetry of the tiling is 5.")

find_rotational_symmetry()