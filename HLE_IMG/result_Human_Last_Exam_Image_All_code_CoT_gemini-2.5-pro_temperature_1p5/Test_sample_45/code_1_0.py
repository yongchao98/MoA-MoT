def find_rotational_symmetry():
    """
    This function determines and explains the rotational symmetry of the given tiling.
    
    Through visual inspection, we can identify several points of rotational symmetry:
    1. The center of each dark blue hexagon is a point of 6-fold symmetry.
    2. The center of each light yellow square is a point of 4-fold symmetry.
    3. The point where the acute angles of two orange rhombuses meet is a point of 2-fold symmetry.

    The overall rotational symmetry of a pattern is defined by the highest order found.
    In this case, the highest order is 6.
    """
    
    # The highest order of rotational symmetry found in the tiling.
    n_fold_symmetry = 6
    
    # Total degrees in a full rotation.
    full_rotation_degrees = 360
    
    # Calculate the minimum angle of rotation that leaves the pattern invariant.
    minimum_angle = full_rotation_degrees / n_fold_symmetry
    
    print(f"The tiling has multiple centers of rotational symmetry.")
    print(f"The highest order of rotational symmetry is {n_fold_symmetry}-fold.")
    print("This corresponds to rotating the tiling around the center of any of the hexagons.")
    print("The calculation for the minimum angle of rotation is:")
    # Print the equation with each number.
    print(f"{full_rotation_degrees} / {n_fold_symmetry} = {int(minimum_angle)} degrees.")

find_rotational_symmetry()

print("<<<6>>>")