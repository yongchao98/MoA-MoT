def find_rotational_symmetry():
    """
    This function explains the process of finding the rotational symmetry
    of the provided tiling through observation.
    """
    
    # Step 1: Identify a center of rotation.
    # We choose the center of one of the prominent dark blue polygons.
    print("Step 1: Identify a potential center of rotation in the tiling.")
    print("A good candidate is the center of a dark blue, seven-sided polygon (heptagon).\n")

    # Step 2: Observe the arrangement of shapes around this center.
    # The central blue heptagon is surrounded by 7 yellow pentagons.
    print("Step 2: Observe the shapes arranged around this central point.")
    print("The central heptagon is surrounded by a ring of 7 identical yellow pentagons.\n")

    # Step 3: Count the number of repeating units to determine the order of symmetry.
    # The number of identical, symmetrically placed units around the center gives us 'n'.
    n = 7
    print(f"Step 3: Count the number of repeating units in this ring.")
    print(f"There are {n} such units. This indicates a {n}-fold rotational symmetry.\n")

    # Step 4: State the conclusion and the corresponding angle of rotation.
    angle = 360 / n
    print("Step 4: Conclude the rotational symmetry and calculate the angle.")
    print(f"The tiling has a {n}-fold rotational symmetry.")
    print("This means rotating the pattern by 360 / 7 degrees results in an identical image.")
    
    # Print the final calculation as requested.
    print("\nFinal calculation:")
    print(f"360 / {n} = {angle:.2f} degrees")

find_rotational_symmetry()