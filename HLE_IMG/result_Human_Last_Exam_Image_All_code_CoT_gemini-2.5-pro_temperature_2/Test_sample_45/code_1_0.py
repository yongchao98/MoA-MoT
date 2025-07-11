def find_rotational_symmetry():
    """
    This function explains and calculates the rotational symmetry of the given tiling.
    """
    print("To find the rotational symmetry of the tiling, we must identify a point of rotation and the smallest angle of rotation that leaves the pattern unchanged.")
    
    # Step 1: Identify a center of rotation
    print("\nLet's choose the center of one of the dark blue, regular hexagons as the point of rotation.")

    # Step 2: Determine the smallest angle of rotation
    print("A regular hexagon has 6 equal sides and angles. By observing the pattern, we can see that if we rotate it around the hexagon's center, the entire tiling matches its original state after a rotation of 1/6th of a full circle.")
    
    full_rotation = 360
    divisions = 6
    smallest_angle = full_rotation / divisions
    
    print(f"The smallest angle of rotation is {full_rotation} / {divisions} = {int(smallest_angle)} degrees.")

    # Step 3: Calculate the order of rotational symmetry
    print("\nThe order of rotational symmetry is 360 degrees divided by the smallest angle of rotation.")
    
    order_of_symmetry = full_rotation / smallest_angle
    
    print(f"\nFinal calculation: {full_rotation} / {int(smallest_angle)} = {int(order_of_symmetry)}")
    print(f"Thus, the tiling has a {int(order_of_symmetry)}-fold rotational symmetry.")

find_rotational_symmetry()