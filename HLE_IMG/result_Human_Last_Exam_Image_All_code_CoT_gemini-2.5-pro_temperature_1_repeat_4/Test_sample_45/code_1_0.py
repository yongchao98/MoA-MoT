def solve_rotational_symmetry():
    """
    This function explains the steps to find the rotational symmetry of the tiling.
    """
    print("To find the rotational symmetry, we look for a center of rotation and count how many times the pattern looks the same during a 360-degree turn.")
    
    print("\nStep 1: Identify a center of rotation.")
    print("Let's choose the center of one of the dark blue, regular hexagons as our point of rotation.")

    print("\nStep 2: Analyze the surrounding pattern.")
    print("Each blue hexagon is surrounded by 6 identical light-yellow pentagonal shapes.")
    
    n_sides = 6
    print(f"The number of identical repeating units around the center is {n_sides}.")

    print("\nStep 3: Calculate the smallest angle of rotation.")
    print("The smallest angle that leaves the pattern unchanged is found by dividing 360 degrees by the number of repeating units.")
    
    angle = 360 / n_sides
    
    print(f"The calculation is: 360 / {n_sides} = {int(angle)} degrees.")

    print("\nStep 4: Determine the order of rotational symmetry.")
    print(f"Since the pattern matches itself {n_sides} times in a full 360-degree rotation, the tiling has a {n_sides}-fold rotational symmetry.")

solve_rotational_symmetry()